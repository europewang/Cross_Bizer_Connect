# -*- coding: utf-8 -*-
import arcpy
import math
import os

def get_distance(p1, p2):
    """计算两个 arcpy.Point 之间的欧氏距离"""
    return math.sqrt((p1.X - p2.X)**2 + (p1.Y - p2.Y)** 2)

def get_line_angle(line):
    """计算线段的角度（弧度）"""
    start_point = line.firstPoint
    end_point = line.lastPoint
    dx = end_point.X - start_point.X
    dy = end_point.Y - start_point.Y
    return math.atan2(dy, dx)

def angle_difference(angle1, angle2):
    """计算两个角度之间的最小差值（弧度）"""
    diff = abs(angle1 - angle2)
    return min(diff, 2 * math.pi - diff)

def classify_lines_by_slope(lines, tolerance_degrees=30):
    """按斜率将线要素分为两组"""
    try:
        if not lines or len(lines) < 2:
            return [], []
            
        # 计算所有线的斜率角度
        line_angles = []
        for i, line in enumerate(lines):
            try:
                start_point = line.firstPoint
                end_point = line.lastPoint
                
                if not start_point or not end_point:
                    continue
                    
                dx = end_point.X - start_point.X
                dy = end_point.Y - start_point.Y
                
                if dx == 0 and dy == 0:
                    continue
                    
                # 计算角度（弧度转角度）
                angle = math.degrees(math.atan2(dy, dx))
                # 标准化到0-180度范围
                if angle < 0:
                    angle += 180
                    
                line_angles.append((i, angle, line))
                
            except Exception as e:
                arcpy.AddWarning(u"计算线{0}斜率失败: {1}".format(i+1, e))
                continue
        
        if len(line_angles) < 2:
            arcpy.AddWarning(u"有效线要素不足，无法分组")
            return [], []
        
        # 使用第一条线作为参考
        reference_angle = line_angles[0][1]
        group_a = [line_angles[0][2]]  # 包含参考线
        group_b = []
        
        tolerance_rad = math.radians(tolerance_degrees)
        
        for i, angle, line in line_angles[1:]:
            # 计算角度差
            angle_diff = abs(angle - reference_angle)
            # 处理角度跨越180度的情况
            if angle_diff > 90:
                angle_diff = 180 - angle_diff
                
            angle_diff_rad = math.radians(angle_diff)
            
            if angle_diff_rad <= tolerance_rad:
                group_a.append(line)
            else:
                group_b.append(line)
        
        arcpy.AddMessage(u"斜率分组完成：参考角度{0:.1f}度，容差{1}度".format(reference_angle, tolerance_degrees))
        return group_a, group_b
        
    except Exception as e:
        arcpy.AddError(u"classify_lines_by_slope函数执行失败: {0}".format(e))
        return [], []

def extend_line(line, length=1000):
    """延长线要素"""
    try:
        if not line:
            return None
            
        start_point = line.firstPoint
        end_point = line.lastPoint
        
        if not start_point or not end_point:
            return None
            
        # 计算方向向量
        dx = end_point.X - start_point.X
        dy = end_point.Y - start_point.Y
        line_length = math.sqrt(dx * dx + dy * dy)
        
        if line_length == 0:
            return None
            
        # 归一化方向向量
        unit_dx = dx / line_length
        unit_dy = dy / line_length
        
        # 向两端延长
        new_start_x = start_point.X - unit_dx * length
        new_start_y = start_point.Y - unit_dy * length
        new_end_x = end_point.X + unit_dx * length
        new_end_y = end_point.Y + unit_dy * length
        
        # 创建延长后的线
        new_start = arcpy.Point(new_start_x, new_start_y)
        new_end = arcpy.Point(new_end_x, new_end_y)
        
        point_array = arcpy.Array([new_start, new_end])
        extended_line = arcpy.Polyline(point_array, line.spatialReference)
        
        return extended_line
        
    except Exception as e:
        arcpy.AddWarning(u"延长线失败: {0}".format(e))
        return None

def find_line_intersection(line1, line2):
    """计算两条线的交点"""
    try:
        if not line1 or not line2:
            return None
            
        # 获取线的端点
        x1, y1 = line1.firstPoint.X, line1.firstPoint.Y
        x2, y2 = line1.lastPoint.X, line1.lastPoint.Y
        x3, y3 = line2.firstPoint.X, line2.firstPoint.Y
        x4, y4 = line2.lastPoint.X, line2.lastPoint.Y
        
        # 计算交点
        denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
        
        if abs(denom) < 1e-10:  # 线平行
            return None
            
        t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom
        
        # 计算交点坐标
        intersection_x = x1 + t * (x2 - x1)
        intersection_y = y1 + t * (y2 - y1)
        
        return arcpy.Point(intersection_x, intersection_y)
        
    except Exception as e:
        arcpy.AddWarning(u"计算交点失败: {0}".format(e))
        return None

def create_bezier_from_endpoints(P0, P3, tan_p0_source, tan_p3_source, fullness, num_points, spatial_ref):
    """
    根据给定的起止点和切线方向点，创建贝塞尔曲线。
    """
    # 限制饱满度因子，防止坐标超出边界
    max_fullness = 0.8  # 限制最大饱满度
    fullness = min(abs(fullness), max_fullness)
    
    # 控制点1，位于起点切线上
    P1_x = P0.X + fullness * (P0.X - tan_p0_source.X)
    P1_y = P0.Y + fullness * (P0.Y - tan_p0_source.Y)
    P1 = arcpy.Point(P1_x, P1_y)

    # 控制点2，位于终点切线上
    P2_x = P3.X + fullness * (P3.X - tan_p3_source.X)
    P2_y = P3.Y + fullness * (P3.Y - tan_p3_source.Y)
    P2 = arcpy.Point(P2_x, P2_y)

    # 检查坐标是否在合理范围内
    def validate_coordinate(coord):
        max_coord = 1e8   # 安全范围
        min_coord = -1e8
        if abs(coord) > max_coord:
            arcpy.AddWarning(u"坐标值 {0} 超出安全范围，已限制到 {1}".format(coord, max_coord if coord > 0 else min_coord))
        return max(min_coord, min(max_coord, coord))
    
    # 生成贝塞尔曲线上的点
    curve_points = []
    for i in range(num_points + 1):
        t = float(i) / num_points
        inv_t = 1 - t
        
        x = (inv_t**3 * P0.X) + (3 * inv_t**2 * t * P1.X) + (3 * inv_t * t**2 * P2.X) + (t**3 * P3.X)
        y = (inv_t**3 * P0.Y) + (3 * inv_t**2 * t * P1.Y) + (3 * inv_t * t**2 * P2.Y) + (t**3 * P3.Y)
        
        x = validate_coordinate(x)
        y = validate_coordinate(y)
        
        curve_points.append(arcpy.Point(x, y))
        
    # 创建折线
    try:
        return arcpy.Polyline(arcpy.Array(curve_points), spatial_ref)
    except Exception as e:
        arcpy.AddWarning(u"创建贝塞尔曲线时出错: {0}".format(e))
        return arcpy.Polyline(arcpy.Array([P0, P3]), spatial_ref)

def get_line_direction_vector(line):
    """
    获取线要素的方向向量
    """
    try:
        if not line or line.length == 0:
            return None
            
        start_point = line.firstPoint
        end_point = line.lastPoint
        
        if not start_point or not end_point:
            return None
            
        # 计算方向向量
        dx = end_point.X - start_point.X
        dy = end_point.Y - start_point.Y
        
        # 归一化
        length = math.sqrt(dx * dx + dy * dy)
        if length == 0:
            return None
            
        return (dx / length, dy / length)
        
    except Exception as e:
        arcpy.AddWarning(u"获取线方向向量失败: {0}".format(e))
        return None

def arrange_lines_by_direction(lines):
    """
    按照线的方向排列线要素，确保能够头连尾连接
    返回: (arranged_lines, failed_lines) - 成功排列的线和无法连接的线
    """
    try:
        if len(lines) < 2:
            return lines, []
        
        arranged = [lines[0]]  # 从第一条线开始
        remaining = lines[1:]
        failed_lines = []
        
        arcpy.AddMessage(u"开始按方向排列线要素，总共 {0} 条线".format(len(lines)))
        
        while remaining:
            current_line = arranged[-1]
            
            # 安全地获取当前线的终点
            try:
                current_end = current_line.lastPoint
                if not current_end:
                    arcpy.AddWarning(u"当前线的终点为空，跳过剩余连接")
                    failed_lines.extend(remaining)
                    break
            except Exception as e:
                arcpy.AddWarning(u"获取当前线终点失败: {0}".format(e))
                failed_lines.extend(remaining)
                break
            
            # 寻找能够连接到当前线尾部的线
            best_match = None
            best_distance = float('inf')
            best_index = -1
            best_reversed = False
            
            for i, candidate_line in enumerate(remaining):
                try:
                    candidate_start = candidate_line.firstPoint
                    candidate_end = candidate_line.lastPoint
                    
                    if not candidate_start or not candidate_end:
                        arcpy.AddWarning(u"候选线 {0} 的端点为空，跳过".format(i))
                        continue
                    
                    # 检查候选线的起点是否接近当前线的终点（正向连接）
                    dist_to_start = get_distance(current_end, candidate_start)
                    # 检查候选线的终点是否接近当前线的终点（反向连接）
                    dist_to_end = get_distance(current_end, candidate_end)
                    
                    # 选择距离更近的连接方式
                    if dist_to_start < best_distance:
                        best_match = candidate_line
                        best_distance = dist_to_start
                        best_index = i
                        best_reversed = False
                    
                    if dist_to_end < best_distance:
                        best_match = candidate_line
                        best_distance = dist_to_end
                        best_index = i
                        best_reversed = True
                        
                except Exception as e:
                    arcpy.AddWarning(u"处理候选线 {0} 时出错: {1}".format(i, e))
                    continue
            
            # 如果找到了合适的连接线
            if best_match and best_distance < 1000:  # 设置一个合理的距离阈值
                try:
                    if best_reversed:
                        arcpy.AddMessage(u"线要素需要反向以实现头尾连接，距离: {0:.2f}".format(best_distance))
                    else:
                        arcpy.AddMessage(u"线要素正向连接，距离: {0:.2f}".format(best_distance))
                    
                    arranged.append(best_match)
                    remaining.pop(best_index)
                    arcpy.AddMessage(u"成功连接线要素，剩余 {0} 条线".format(len(remaining)))
                except Exception as e:
                    arcpy.AddWarning(u"添加连接线时出错: {0}".format(e))
                    failed_lines.extend(remaining)
                    break
            else:
                # 无法找到合适的连接线，将剩余的线标记为失败
                arcpy.AddWarning(u"无法为当前线找到合适的头尾连接（最近距离: {0:.2f}），剩余 {1} 条线无法连接".format(best_distance, len(remaining)))
                failed_lines.extend(remaining)
                break
        
        arcpy.AddMessage(u"线要素排列完成：成功排列 {0} 条，失败 {1} 条".format(len(arranged), len(failed_lines)))
        return arranged, failed_lines
        
    except Exception as e:
        arcpy.AddError(u"arrange_lines_by_direction函数执行失败: {0}".format(e))
        import traceback
        arcpy.AddError(traceback.format_exc())
        return lines[:1] if lines else [], lines[1:] if len(lines) > 1 else []

def create_directional_bezier_connections(lines, fullness_factor, num_points, spatial_ref):
    """
    按照线的方向创建头尾连接的贝塞尔曲线
    """
    try:
        if len(lines) < 2:
            arcpy.AddWarning(u"线要素数量不足，无法创建连接")
            return []
        
        curves = []
        arcpy.AddMessage(u"开始创建 {0} 条线之间的方向性贝塞尔曲线连接".format(len(lines)))
        
        for i in range(len(lines) - 1):
            try:
                current_line = lines[i]
                next_line = lines[i + 1]
                
                # 安全地获取线的端点
                try:
                    current_end = current_line.lastPoint
                    current_start = current_line.firstPoint
                    next_start = next_line.firstPoint
                    next_end = next_line.lastPoint
                    
                    if not all([current_end, current_start, next_start, next_end]):
                        arcpy.AddWarning(u"线{0}或线{1}的端点为空，跳过连接".format(i+1, i+2))
                        continue
                        
                except Exception as e:
                    arcpy.AddWarning(u"获取线{0}或线{1}的端点失败: {2}".format(i+1, i+2, e))
                    continue
                
                # 检查连接距离
                try:
                    connection_distance = get_distance(current_end, next_start)
                    reverse_distance = get_distance(current_end, next_end)
                    
                    arcpy.AddMessage(u"线{0}到线{1}: 正向距离={2:.2f}, 反向距离={3:.2f}".format(
                        i+1, i+2, connection_distance, reverse_distance))
                        
                except Exception as e:
                    arcpy.AddWarning(u"计算线{0}到线{1}的连接距离失败: {2}".format(i+1, i+2, e))
                    continue
                
                if reverse_distance < connection_distance:
                    # 使用反向连接（当前线尾连接下一条线尾）
                    start_point = current_end
                    end_point = next_end
                    # 计算切线方向
                    current_tangent = current_start  # 当前线的起点作为切线参考
                    next_tangent = next_start  # 下一条线的起点作为切线参考
                    arcpy.AddMessage(u"使用反向连接: 线{0}尾 -> 线{1}尾，距离: {2:.2f}".format(i+1, i+2, reverse_distance))
                else:
                    # 使用正向连接（当前线尾连接下一条线头）
                    start_point = current_end
                    end_point = next_start
                    # 计算切线方向
                    current_tangent = current_start  # 当前线的起点作为切线参考
                    next_tangent = next_end  # 下一条线的终点作为切线参考
                    arcpy.AddMessage(u"使用正向连接: 线{0}尾 -> 线{1}头，距离: {2:.2f}".format(i+1, i+2, connection_distance))
                
                # 创建贝塞尔曲线
                try:
                    curve = create_bezier_from_endpoints(
                        start_point, end_point, 
                        current_tangent, next_tangent,
                        fullness_factor, num_points, spatial_ref
                    )
                    if curve:
                        curves.append(curve)
                        arcpy.AddMessage(u"成功创建线{0}到线{1}的方向性贝塞尔曲线".format(i+1, i+2))
                    else:
                        arcpy.AddWarning(u"线{0}到线{1}的贝塞尔曲线创建返回空值".format(i+1, i+2))
                except Exception as e:
                    arcpy.AddWarning(u"创建线{0}到线{1}的贝塞尔曲线失败: {2}".format(i+1, i+2, e))
                    import traceback
                    arcpy.AddWarning(traceback.format_exc())
                    
            except Exception as e:
                arcpy.AddWarning(u"处理线{0}到线{1}的连接时出错: {2}".format(i+1, i+2, e))
                continue
        
        arcpy.AddMessage(u"方向性贝塞尔曲线创建完成，成功创建 {0} 条曲线".format(len(curves)))
        return curves
        
    except Exception as e:
        arcpy.AddError(u"create_directional_bezier_connections函数执行失败: {0}".format(e))
        import traceback
        arcpy.AddError(traceback.format_exc())
        return []

def create_cross_group_bezier_connections(group_a, group_b, fullness_factor, num_points, spatial_ref):
    """
    在A组和B组线要素之间创建贝塞尔曲线连接
    找到A线和B线所有端点中最近的两个点，根据端点类型确定曲线方向
    """
    try:
        curves = []
        arcpy.AddMessage(u"开始创建A组({0}条)和B组({1}条)之间的贝塞尔曲线连接".format(len(group_a), len(group_b)))
        
        processed_pairs = set()
        
        # 为每条A线找到最近的B线端点
        for i, line_a in enumerate(group_a):
            for j, line_b in enumerate(group_b):
                # 跳过已处理的配对
                if (i, j) in processed_pairs:
                    continue
                
                # 收集A线和B线的所有端点
                a_points = [
                    (line_a.firstPoint, 'a_start', line_a),
                    (line_a.lastPoint, 'a_end', line_a)
                ]
                b_points = [
                    (line_b.firstPoint, 'b_start', line_b),
                    (line_b.lastPoint, 'b_end', line_b)
                ]
                
                # 找到所有端点中距离最近的两个点
                min_distance = float('inf')
                best_connection = None
                
                for a_point, a_type, a_line in a_points:
                    for b_point, b_type, b_line in b_points:
                        distance = get_distance(a_point, b_point)
                        if distance < min_distance:
                            min_distance = distance
                            best_connection = (a_point, b_point, a_type, b_type, a_line, b_line)
                
                if best_connection:
                    a_point, b_point, a_type, b_type, a_line, b_line = best_connection
                    
                    # 标记这对线已处理
                    processed_pairs.add((i, j))
                    
                    # 根据端点类型确定贝塞尔曲线方向和是否符合头尾相连原则
                    curve_start = None
                    curve_end = None
                    connection_valid = False
                    
                    if a_type == 'a_start' and b_type == 'b_end':
                        # A线起点连B线终点：符合头尾相连原则
                        curve_start = a_point
                        curve_end = b_point
                        connection_valid = True
                        arcpy.AddMessage(u"A组线{0}起点连接B组线{1}终点".format(i+1, j+1))
                    elif a_type == 'a_end' and b_type == 'b_start':
                        # A线终点连B线起点：符合头尾相连原则
                        curve_start = b_point
                        curve_end = a_point
                        connection_valid = True
                        arcpy.AddMessage(u"B组线{0}起点连接A组线{1}终点".format(j+1, i+1))
                    elif a_type == 'a_start' and b_type == 'b_start':
                        # A线起点连B线起点：不符合头尾相连原则
                        arcpy.AddWarning(u"A组线{0}起点连接B组线{1}起点，不符合头尾相连原则，跳过".format(i+1, j+1))
                        continue
                    elif a_type == 'a_end' and b_type == 'b_end':
                        # A线终点连B线终点：不符合头尾相连原则
                        arcpy.AddWarning(u"A组线{0}终点连接B组线{1}终点，不符合头尾相连原则，跳过".format(i+1, j+1))
                        continue
                    
                    if connection_valid and curve_start and curve_end:
                        # 延长两条线并找到交点作为控制点
                        try:
                            extended_a = extend_line(a_line, length=1000)
                            extended_b = extend_line(b_line, length=1000)
                            
                            if extended_a and extended_b:
                                intersection = find_line_intersection(extended_a, extended_b)
                                
                                if intersection:
                                    # 使用交点作为控制点创建贝塞尔曲线
                                    curve = create_bezier_curve_with_control_point(
                                        curve_start, curve_end, intersection, 
                                        fullness_factor, num_points, spatial_ref
                                    )
                                    
                                    if curve:
                                        curves.append(curve)
                                        arcpy.AddMessage(u"成功创建贝塞尔曲线，使用延长线交点作为控制点")
                                    else:
                                        arcpy.AddWarning(u"贝塞尔曲线创建失败")
                                else:
                                    # 如果没有交点，使用中点作为控制点
                                    mid_x = (curve_start.X + curve_end.X) / 2
                                    mid_y = (curve_start.Y + curve_end.Y) / 2
                                    mid_point = arcpy.Point(mid_x, mid_y)
                                    
                                    curve = create_bezier_curve_with_control_point(
                                        curve_start, curve_end, mid_point,
                                        fullness_factor, num_points, spatial_ref
                                    )
                                    
                                    if curve:
                                        curves.append(curve)
                                        arcpy.AddMessage(u"使用中点作为控制点创建贝塞尔曲线")
                            else:
                                arcpy.AddWarning(u"线延长失败")
                                
                        except Exception as curve_error:
                            arcpy.AddWarning(u"创建贝塞尔曲线时出错: {0}".format(curve_error))
                            continue
        
        arcpy.AddMessage(u"AB组贝塞尔曲线连接完成，成功创建 {0} 条曲线".format(len(curves)))
        return curves
        
    except Exception as e:
        arcpy.AddError(u"create_cross_group_bezier_connections函数执行失败: {0}".format(e))
        return []

def create_bezier_curve_with_control_point(start_point, end_point, control_point, fullness_factor, num_points, spatial_ref):
    """
    使用指定的控制点创建贝塞尔曲线
    """
    try:
        # 创建贝塞尔曲线点
        bezier_points = []
        for i in range(num_points + 1):
            t = float(i) / num_points
            
            # 二次贝塞尔曲线公式 (使用单个控制点)
            x = ((1-t)**2 * end_point.X + 
                 2*(1-t)*t * control_point.X + 
                 t**2 * start_point.X)
            
            y = ((1-t)**2 * end_point.Y + 
                 2*(1-t)*t * control_point.Y + 
                 t**2 * start_point.Y)
            
            bezier_points.append(arcpy.Point(x, y))
        
        # 创建折线几何
        polyline = arcpy.Polyline(arcpy.Array(bezier_points), spatial_ref)
        return polyline
        
    except Exception as e:
        arcpy.AddWarning(u"创建贝塞尔曲线时出错: {0}".format(e))
        return None

if __name__ == '__main__':
    try:
        arcpy.AddMessage(u"脚本开始执行...")

        # 获取参数
        arcpy.AddMessage(u"正在获取参数 0 (输入要素)...")
        input_features = arcpy.GetParameter(0)
        if not input_features:
            raise ValueError(u"输入要素 (参数 0) 未提供或无效。")
        arcpy.AddMessage(u"参数 0 获取成功。")

        # 不再需要输出路径参数，直接写入输入要素类
        output_fc = input_features
        arcpy.AddMessage(u"将在输入要素类中添加贝塞尔曲线")

        arcpy.AddMessage(u"正在获取参数 1 (曲线饱满度)...")
        fullness_factor = arcpy.GetParameter(1)
        if fullness_factor is None:
            raise ValueError(u"曲线饱满度 (参数 1) 未提供。")
        arcpy.AddMessage(u"参数 1 (曲线饱满度): {}".format(fullness_factor))

        arcpy.AddMessage(u"正在获取参数 2 (曲线平滑度)...")
        num_points = arcpy.GetParameter(2)
        if num_points is None or num_points < 2:
            num_points = 20  # 默认值
            arcpy.AddWarning(u"曲线平滑度参数无效，使用默认值: {}".format(num_points))
        arcpy.AddMessage(u"参数 2 (曲线平滑度): {}".format(num_points))

        arcpy.AddMessage(u"所有参数获取成功。")

        # 输入要素验证
        arcpy.AddMessage(u"正在验证输入要素...")
        desc = arcpy.Describe(input_features)
        if desc.shapeType != "Polyline":
            raise ValueError(u"输入要素必须是线要素图层（Polyline），当前类型为：{0}".format(desc.shapeType))
        arcpy.AddMessage(u"输入要素类型验证通过: 线要素图层")

        # 确定要素ID字段名（解决Shapefile和GDB的差异）# 确定要素ID字段
        id_field = "FID"
        # 对于不同数据源，可能需要使用不同的ID字段
        if hasattr(desc, 'catalogPath') and desc.catalogPath and desc.catalogPath.endswith('.shp'):
            id_field = "FID"  # Shapefile使用FID
        elif desc.dataType == "FeatureClass":
            id_field = "OBJECTID"  # 地理数据库要素类使用OBJECTID
        arcpy.AddMessage(u"使用的要素ID字段: {0}".format(id_field))

        # 分析输入要素信息
        arcpy.AddMessage(u"正在分析输入要素...")
        arcpy.AddMessage(u"输入要素类型: {0}".format(desc.dataType))
        arcpy.AddMessage(u"输入要素路径: {0}".format(desc.catalogPath if hasattr(desc, 'catalogPath') else '未知'))
        
        # 获取空间参考
        spatial_ref = desc.spatialReference
        arcpy.AddMessage(u"空间参考: {0}".format(spatial_ref.name if spatial_ref else '未定义'))
        
        # # 处理未定义的空间参考
        # if not spatial_ref or spatial_ref.name == "Unknown":
        #     arcpy.AddWarning(u"检测到未定义的空间参考系统，将使用默认的地理坐标系统")
        #     try:
        #         spatial_ref = arcpy.SpatialReference(4326)  # WGS84
        #         arcpy.AddMessage(u"已设置空间参考为: {0}".format(spatial_ref.name))
        #     except:
        #         try:
        #             spatial_ref = arcpy.SpatialReference(3857)  # Web Mercator
        #             arcpy.AddMessage(u"已设置空间参考为: {0}".format(spatial_ref.name))
        #         except:
        #             spatial_ref = None
        #             arcpy.AddWarning(u"无法设置空间参考系统，将使用无空间参考模式")

        # 直接使用输入要素类，不需要创建新的输出要素类
        arcpy.AddMessage(u"将直接在输入要素类中添加贝塞尔曲线")

        # 读取输入的线要素
        arcpy.AddMessage(u"正在读取输入要素...")
        lines = []
        feature_count = 0
        total_count = 0
        
        # 获取要素总数
        try:
            result = arcpy.GetCount_management(input_features)
            total_count = int(result.getOutput(0))
            arcpy.AddMessage(u"输入要素总数: {0}".format(total_count))
        except:
            arcpy.AddMessage(u"无法获取要素总数，继续处理...")
        
        # 检查选择集
        has_user_selection = False
        selection_oids = []
        if hasattr(desc, 'FIDSet') and desc.FIDSet:
            arcpy.AddMessage(u"有选择集：FIDSet: {0}".format(desc.FIDSet))
            selection_oids = [int(fid) for fid in desc.FIDSet.split(';') if fid.strip()]
            if len(selection_oids) > 0:
                has_user_selection = True
                arcpy.AddMessage(u"检测到选择集，包含 {0} 个要素".format(len(selection_oids)))

        # 读取要素 - 优先处理选择集
        if has_user_selection and len(selection_oids) > 0:
            arcpy.AddMessage(u"【debug1】使用CopyFeatures方法处理选中要素...")
            try:
                # 创建临时要素类来存储选择集，使用系统临时目录
                import tempfile
                temp_dir = tempfile.gettempdir()
                temp_fc = os.path.join(temp_dir, "temp_selected_features.shp")
                arcpy.AddMessage(u"【debug】创建临时要素类: {0}".format(temp_fc))
                
                # 删除可能存在的临时文件
                if arcpy.Exists(temp_fc):
                    arcpy.Delete_management(temp_fc)
                
                # 使用CopyFeatures复制选择集到临时要素类
                arcpy.CopyFeatures_management(input_features, temp_fc)
                arcpy.AddMessage(u"【debug】CopyFeatures执行成功")
                
                # 从临时要素类读取要素
                arcpy.AddMessage(u"【debug】从临时要素类读取要素...")
                with arcpy.da.SearchCursor(temp_fc, ["SHAPE@", "OID@"]) as cursor:
                    arcpy.AddMessage(u"【debug】SearchCursor创建成功，开始遍历...")
                    for row in cursor:
                        geometry = row[0]
                        oid = row[1]
                        feature_count += 1
                        arcpy.AddMessage(u"【debug】读取到要素 OID: {0}".format(oid))
                        
                        if geometry is not None:
                            try:
                                if geometry.type.lower() == "polyline":
                                    if geometry.length > 0 and geometry.pointCount > 1:
                                        lines.append(geometry)
                                        arcpy.AddMessage(u"读取到有效线要素 OID: {0}, 长度: {1:.2f}".format(oid, geometry.length))
                                    else:
                                        arcpy.AddWarning(u"跳过无效线要素 OID: {0} (长度为0或点数不足)".format(oid))
                                else:
                                    arcpy.AddWarning(u"跳过非线要素 OID: {0}, 类型: {1}".format(oid, geometry.type))
                            except Exception as geom_error:
                                arcpy.AddWarning(u"处理几何要素 OID {0} 时出错: {1}".format(oid, geom_error))
                        else:
                            arcpy.AddWarning(u"跳过空几何要素 OID: {0}".format(oid))
                
                arcpy.AddMessage(u"从选择集读取了 {0} 个要素".format(len(lines)))
                
                # 清理临时要素类
                try:
                    if arcpy.Exists(temp_fc):
                        arcpy.Delete_management(temp_fc)
                        arcpy.AddMessage(u"【debug】临时要素类已清理")
                except:
                    pass
                    
            except Exception as cursor_error:
                arcpy.AddWarning(u"处理选择集失败: {0}".format(cursor_error))
                arcpy.AddMessage(u"回退到处理所有要素...")
                
                # 清理临时要素类
                try:
                    if arcpy.Exists(temp_fc):
                        arcpy.Delete_management(temp_fc)
                except:
                    pass
                    
                has_user_selection = False  # 失败则处理所有要素
                lines = []  # 清空已读取的数据
                feature_count = 0

        # 如果没有选择集或选择集处理失败，读取所有要素
        if not has_user_selection or len(lines) == 0:
            arcpy.AddMessage(u"读取所有线要素...")
            try:
                # 优先使用catalogPath，如果不存在则使用input_features
                data_source = desc.catalogPath if hasattr(desc, 'catalogPath') and desc.catalogPath else input_features
                arcpy.AddMessage(u"【debug】fallback使用数据源: {0}".format(data_source))
                arcpy.AddMessage(u"【debug】尝试使用SearchCursor读取所有要素...")
                with arcpy.da.SearchCursor(data_source, ["SHAPE@", "OID@"]) as cursor:
                    arcpy.AddMessage(u"【debug】SearchCursor创建成功，开始遍历...")
                    for row in cursor:
                         geometry = row[0]
                         oid = row[1]
                         feature_count += 1
                         arcpy.AddMessage(u"【debug】读取到要素 OID: {0}".format(oid))
                         
                         if geometry is not None:
                             try:
                                 if geometry.type.lower() == "polyline":
                                     if geometry.length > 0 and geometry.pointCount > 1:
                                         lines.append(geometry)
                                         arcpy.AddMessage(u"读取到有效线要素 OID: {0}, 长度: {1:.2f}".format(oid, geometry.length))
                                     else:
                                         arcpy.AddWarning(u"跳过无效线要素 OID: {0} (长度为0或点数不足)".format(oid))
                                 else:
                                     arcpy.AddWarning(u"跳过非线要素 OID: {0}, 类型: {1}".format(oid, geometry.type))
                             except Exception as geom_error:
                                 arcpy.AddWarning(u"处理几何要素 OID {0} 时出错: {1}".format(oid, geom_error))
                         else:
                             arcpy.AddWarning(u"跳过空几何要素 OID: {0}".format(oid))
            except Exception as cursor_error:
                arcpy.AddWarning(u"读取要素失败: {0}".format(cursor_error))
                # 如果还是失败，尝试使用另一种数据源
                try:
                    alternative_source = input_features if data_source != input_features else desc.catalogPath
                    if alternative_source:
                        arcpy.AddMessage(u"【debug】尝试使用备用数据源: {0}".format(alternative_source))
                        with arcpy.da.SearchCursor(alternative_source, ["SHAPE@", "OID@"]) as cursor:
                            for row in cursor:
                                geometry = row[0]
                                oid = row[1]
                                feature_count += 1
                                
                                if geometry is not None:
                                    try:
                                        if geometry.type.lower() == "polyline":
                                            if geometry.length > 0 and geometry.pointCount > 1:
                                                lines.append(geometry)
                                                arcpy.AddMessage(u"读取到有效线要素 OID: {0}, 长度: {1:.2f}".format(oid, geometry.length))
                                            else:
                                                arcpy.AddWarning(u"跳过无效线要素 OID: {0} (长度为0或点数不足)".format(oid))
                                        else:
                                            arcpy.AddWarning(u"跳过非线要素 OID: {0}, 类型: {1}".format(oid, geometry.type))
                                    except Exception as geom_error:
                                        arcpy.AddWarning(u"处理几何要素 OID {0} 时出错: {1}".format(oid, geom_error))
                                else:
                                    arcpy.AddWarning(u"跳过空几何要素 OID: {0}".format(oid))
                except Exception as alternative_error:
                    arcpy.AddWarning(u"使用备用数据源读取要素也失败: {0}".format(alternative_error))
        
        arcpy.AddMessage(u"总共处理了 {0} 个要素，其中 {1} 个是有效的线要素。".format(feature_count, len(lines)))

        if not lines:
            if feature_count == 0:
                raise ValueError(u"输入要素集为空。请确保图层包含线要素且有数据。")
            else:
                raise ValueError(u"输入的 {0} 个要素中没有找到有效线要素。".format(feature_count))

        # 使用斜率容差分类线要素为AB两组
        arcpy.AddMessage(u"开始按斜率分类线要素...")
        
        try:
            # 按斜率分类线要素
            group_a, group_b = classify_lines_by_slope(lines, tolerance_degrees=30)
            
            if len(group_a) == 0 or len(group_b) == 0:
                arcpy.AddWarning(u"无法将线要素分为两组，A组: {0} 条，B组: {1} 条".format(len(group_a), len(group_b)))
                arcpy.AddMessage(u"回退到原有的成对连接方式...")
                raise Exception(u"线要素分组失败")
            
            arcpy.AddMessage(u"成功分组：A组 {0} 条线，B组 {1} 条线".format(len(group_a), len(group_b)))
            
            # 创建AB组之间的贝塞尔曲线连接
            curves = create_cross_group_bezier_connections(
                group_a, group_b, fullness_factor, num_points, spatial_ref
            )
            
            if not curves or len(curves) == 0:
                arcpy.AddWarning(u"AB组连接未能创建任何曲线，回退到原有连接方式")
                raise Exception(u"AB组连接创建曲线失败")
                
        except Exception as group_error:
            arcpy.AddWarning(u"AB组连接失败: {0}".format(group_error))
            arcpy.AddMessage(u"回退到原有的成对连接方式...")
            
            # 回退到原有的成对连接方式
            curves = []
            for i in range(len(lines) - 1):
                line1 = lines[i]
                line2 = lines[i + 1]
                
                # 计算四种端点组合的距离和
                distances = [
                    get_distance(line1.firstPoint, line2.firstPoint),
                    get_distance(line1.firstPoint, line2.lastPoint),
                    get_distance(line1.lastPoint, line2.firstPoint),
                    get_distance(line1.lastPoint, line2.lastPoint)
                ]
                
                min_index = distances.index(min(distances))
                
                if min_index == 0:
                    start_point, end_point = line1.firstPoint, line2.firstPoint
                    start_tangent, end_tangent = line1.lastPoint, line2.lastPoint
                elif min_index == 1:
                    start_point, end_point = line1.firstPoint, line2.lastPoint
                    start_tangent, end_tangent = line1.lastPoint, line2.firstPoint
                elif min_index == 2:
                    start_point, end_point = line1.lastPoint, line2.firstPoint
                    start_tangent, end_tangent = line1.firstPoint, line2.lastPoint
                else:
                    start_point, end_point = line1.lastPoint, line2.lastPoint
                    start_tangent, end_tangent = line1.firstPoint, line2.firstPoint
                
                curve = create_bezier_from_endpoints(
                    start_point, end_point, start_tangent, end_tangent,
                    fullness_factor, num_points, spatial_ref
                )
                
                if curve:
                    curves.append(curve)
                    arcpy.AddMessage(u"成功创建第 {0} 条贝塞尔曲线".format(i + 1))
        
        # 将曲线写入输入要素类（使用编辑会话避免锁定问题）
        arcpy.AddMessage(u"开始将贝塞尔曲线写入输入要素类...")
        
        # 获取工作空间路径
        workspace = os.path.dirname(desc.catalogPath) if hasattr(desc, 'catalogPath') and desc.catalogPath else None
        
        if workspace and workspace.endswith('.gdb'):
            # 地理数据库需要编辑会话
            arcpy.AddMessage(u"检测到地理数据库，启动编辑会话...")
            edit = arcpy.da.Editor(workspace)
            edit.startEditing(False, True)
            edit.startOperation()
            
            try:
                with arcpy.da.InsertCursor(input_features, ["SHAPE@"]) as insert_cursor:
                    curves_created = 0
                    for curve in curves:
                        if curve:
                            insert_cursor.insertRow([curve])
                            curves_created += 1
                    
                    arcpy.AddMessage(u"总共成功创建了 {0} 条方向性贝塞尔曲线".format(curves_created))
                
                edit.stopOperation()
                edit.stopEditing(True)
                arcpy.AddMessage(u"编辑会话已保存并关闭")
                
            except Exception as edit_error:
                edit.stopOperation()
                edit.stopEditing(False)
                arcpy.AddError(u"编辑会话中出错: {0}".format(edit_error))
                raise
        else:
             # Shapefile或其他格式，需要特殊处理避免锁定
             arcpy.AddMessage(u"检测到Shapefile格式，使用特殊方法插入...")
             
             # 先清除选择集，避免锁定冲突
             try:
                 arcpy.SelectLayerByAttribute_management(input_features, "CLEAR_SELECTION")
                 arcpy.AddMessage(u"已清除选择集")
             except:
                 arcpy.AddMessage(u"无法清除选择集，继续尝试插入")
             
             # 获取数据源路径
             data_source = desc.catalogPath if hasattr(desc, 'catalogPath') and desc.catalogPath else input_features
             arcpy.AddMessage(u"使用数据源: {0}".format(data_source))
             
             try:
                 with arcpy.da.InsertCursor(data_source, ["SHAPE@"]) as insert_cursor:
                     curves_created = 0
                     for curve in curves:
                         if curve:
                             insert_cursor.insertRow([curve])
                             curves_created += 1
                     
                     arcpy.AddMessage(u"总共成功创建了 {0} 条方向性贝塞尔曲线".format(curves_created))
             except Exception as insert_error:
                 arcpy.AddMessage(u"直接插入失败，尝试使用临时要素类方法: {0}".format(insert_error))
                 
                 # 创建临时要素类存储贝塞尔曲线
                 import tempfile
                 temp_dir = tempfile.gettempdir()
                 temp_curves_fc = os.path.join(temp_dir, "temp_bezier_curves.shp")
                 
                 if arcpy.Exists(temp_curves_fc):
                     arcpy.Delete_management(temp_curves_fc)
                 
                 # 创建临时要素类
                 arcpy.CreateFeatureclass_management(
                     temp_dir, "temp_bezier_curves.shp", "POLYLINE", 
                     spatial_reference=spatial_ref
                 )
                 
                 # 将曲线写入临时要素类
                 with arcpy.da.InsertCursor(temp_curves_fc, ["SHAPE@"]) as temp_cursor:
                     curves_created = 0
                     for curve in curves:
                         if curve:
                             temp_cursor.insertRow([curve])
                             curves_created += 1
                 
                 # 将临时要素类追加到原始要素类
                 arcpy.Append_management(temp_curves_fc, data_source, "NO_TEST")
                 arcpy.AddMessage(u"通过临时要素类成功添加了 {0} 条贝塞尔曲线".format(curves_created))
                 
                 # 清理临时文件
                 try:
                     if arcpy.Exists(temp_curves_fc):
                         arcpy.Delete_management(temp_curves_fc)
                 except:
                     pass

        arcpy.AddMessage(u"贝塞尔曲线创建成功！已添加到输入要素类中")

    except Exception as e:
        arcpy.AddError(u"脚本执行出错: {0}".format(e))
        import traceback
        arcpy.AddError(traceback.format_exc())
