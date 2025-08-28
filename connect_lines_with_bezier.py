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
        
        # 确保A组是线段数量较少的那组
        if len(group_a) > len(group_b):
            arcpy.AddMessage(u"A组线段数({0})>B组线段数({1})，交换A、B组".format(len(group_a), len(group_b)))
            group_a, group_b = group_b, group_a
        else:
            arcpy.AddMessage(u"A组线段数({0})<=B组线段数({1})，保持不变".format(len(group_a), len(group_b)))
        
        return group_a, group_b
        
    except Exception as e:
        arcpy.AddError(u"classify_lines_by_slope函数执行失败: {0}".format(e))
        return [], []

def classify_lines_by_distance(lines):
    """基于线段间距离进行CD分组，用于直线连接
    以第一条线的中点为参考，计算其他线与该线的距离，按距离排序后寻找突变点分组
    """
    try:
        if not lines or len(lines) < 2:
            arcpy.AddWarning(u"线段数量不足，无法进行距离分组")
            return [], []
        
        arcpy.AddMessage(u"开始基于距离进行CD分组...")
        
        # 获取有效线段和中心点
        valid_lines = []
        for line in lines:
            if line and line.centroid:
                valid_lines.append(line)
        
        if len(valid_lines) < 2:
            arcpy.AddWarning(u"有效线段数量不足，无法进行距离分组")
            return [], []
        
        # 以第一条线为参考线
        reference_line = valid_lines[0]
        reference_center = reference_line.centroid
        
        # 计算其他线与参考线的距离
        line_distances = []
        for i, line in enumerate(valid_lines[1:], 1):  # 从第二条线开始
            dist = get_distance(reference_center, line.centroid)
            line_distances.append((dist, i, line))
        
        # 按距离排序
        line_distances.sort()
        
        # 寻找距离突变点
        split_index = find_distance_jump(line_distances)
        
        # 根据突变点分组
        group_c = [reference_line]  # C组包含参考线
        group_d = []
        
        # 将距离较近的线段分到C组，距离较远的分到D组
        for i, (dist, line_idx, line) in enumerate(line_distances):
            if i < split_index:
                group_c.append(line)
            else:
                group_d.append(line)
        
        # 如果没有找到合适的突变点，使用中位数分割
        if not group_d:
            mid_point = len(line_distances) // 2
            group_c = [reference_line]
            group_d = []
            for i, (dist, line_idx, line) in enumerate(line_distances):
                if i < mid_point:
                    group_c.append(line)
                else:
                    group_d.append(line)
        
        # 确保条数多的为D组，条数少的为C组
        if len(group_c) > len(group_d):
            arcpy.AddMessage(u"C组线段数({0})>D组线段数({1})，交换C、D组，确保条数多的为D组".format(len(group_c), len(group_d)))
            group_c, group_d = group_d, group_c
        
        arcpy.AddMessage(u"距离分组完成: C组{0}条线，D组{1}条线".format(len(group_c), len(group_d)))
        
        return group_c, group_d
        
    except Exception as e:
        arcpy.AddError(u"classify_lines_by_distance函数执行失败: {0}".format(e))
        return [], []

def find_distance_jump(line_distances):
    """寻找距离序列中的突变点
    
    Args:
        line_distances: 按距离排序的线段列表 [(distance, index, line), ...]
    
    Returns:
        int: 突变点的索引，用于分组
    """
    if len(line_distances) < 2:
        return 0
    
    # 计算相邻距离的差值
    distance_diffs = []
    for i in range(1, len(line_distances)):
        diff = line_distances[i][0] - line_distances[i-1][0]
        distance_diffs.append(diff)
    
    if not distance_diffs:
        return len(line_distances) // 2
    
    # 寻找最大的距离跳跃
    max_diff = max(distance_diffs)
    max_diff_index = distance_diffs.index(max_diff)
    
    # 只有当最大跳跃显著大于平均跳跃时才认为是突变点
    avg_diff = sum(distance_diffs) / len(distance_diffs)
    
    if max_diff > avg_diff * 2:  # 突变阈值：大于平均值的2倍
        return max_diff_index + 1
    else:
        # 没有明显突变，使用中位数分割
        return len(line_distances) // 2

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

def extend_line_from_last_segment(line, length=1000, use_end_segment=True):
    """基于折线头部或尾部一小段的方向延长线要素
    
    Args:
        line: 输入线要素
        length: 延长长度
        use_end_segment: True=使用尾部一小段，False=使用头部一小段
    """
    try:
        if not line or line.pointCount < 2:
            return None
        
        # 获取折线的所有点
        points = []
        for part in line:
            for point in part:
                points.append(point)
        
        if len(points) < 2:
            return None
        
        if use_end_segment:
            # 使用尾部一小段（最后两个点）
            last_point = points[len(points)-1]  # 最后一个点
            second_last_point = points[len(points)-2]  # 倒数第二个点
            # 打印最后一个点的坐标，保留6位小数
            arcpy.AddMessage(u"最后一个点坐标: X={:.6f}, Y={:.6f}".format(last_point.X, last_point.Y))
            # 打印倒数第二个点的坐标，保留6位小数
            arcpy.AddMessage(u"倒数第二个点坐标: X={:.6f}, Y={:.6f}".format(second_last_point.X, second_last_point.Y))

            # 计算尾部一小段的方向向量
            dx = last_point.X - second_last_point.X
            dy = last_point.Y - second_last_point.Y
            segment_length = math.sqrt(dx * dx + dy * dy)
            
            if segment_length == 0:
                return None
            
            # 归一化方向向量
            unit_dx = dx / segment_length
            unit_dy = dy / segment_length
            
            # 从最后一个点向前延长（基于尾部一小段的方向）
            extended_x = last_point.X + unit_dx * length
            extended_y = last_point.Y + unit_dy * length
            
            # 创建延长后的线（从倒数第二个点到延长点）
            start_point = arcpy.Point(second_last_point.X, second_last_point.Y)
            end_point = arcpy.Point(extended_x, extended_y)
        else:
            # 使用头部一小段（前两个点）
            first_point = points[0]  # 第一个点
            second_point = points[1]  # 第二个点
            # 打印最后一个点的坐标，保留6位小数
            arcpy.AddMessage(u"第一个点坐标: X={:.6f}, Y={:.6f}".format(first_point.X, first_point.Y))
            # 打印倒数第二个点的坐标，保留6位小数
            arcpy.AddMessage(u"第二个点坐标: X={:.6f}, Y={:.6f}".format(second_point.X, second_point.Y))

            
            
            # 计算头部一小段的方向向量
            dx = second_point.X - first_point.X
            dy = second_point.Y - first_point.Y
            segment_length = math.sqrt(dx * dx + dy * dy)
            
            if segment_length == 0:
                return None
            
            # 归一化方向向量
            unit_dx = dx / segment_length
            unit_dy = dy / segment_length
            
            # 从第一个点向后延长（基于头部一小段的反方向）
            extended_x = first_point.X - unit_dx * length
            extended_y = first_point.Y - unit_dy * length
            
            # 创建延长后的线（从延长点到第二个点）
            start_point = arcpy.Point(extended_x, extended_y)
            end_point = arcpy.Point(second_point.X, second_point.Y)
        
        point_array = arcpy.Array([start_point, end_point])
        extended_line = arcpy.Polyline(point_array, line.spatialReference)
        
        return extended_line
        
    except Exception as e:
        arcpy.AddWarning(u"基于线段一小段延长线失败: {0}".format(e))
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

def create_cross_group_bezier_connections(group_a, group_b, fullness_factor, num_points, spatial_ref):
    """
    在A组和B组线要素之间创建贝塞尔曲线连接
    新逻辑：循环找最近点对生成贝塞尔曲线，直到A组剩最后一条线时与B组所有剩余线全连接
    """
    try:
        curves = []
        arcpy.AddMessage(u"开始创建A组({0}条)和B组({1}条)之间的贝塞尔曲线连接".format(len(group_a), len(group_b)))
        
        # 创建可变的线段列表副本
        remaining_a = list(group_a)
        remaining_b = list(group_b)
        
        connection_count = 0
        
        # 循环连接，直到A组剩下最后一条线
        while len(remaining_a) > 1 and len(remaining_b) > 0:
            # 找到A组和B组中距离最近的两个端点
            min_distance = float('inf')
            best_connection = None
            best_a_idx = -1
            best_b_idx = -1
            
            for i, line_a in enumerate(remaining_a):
                for j, line_b in enumerate(remaining_b):
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
                    for a_point, a_type, a_line in a_points:
                        for b_point, b_type, b_line in b_points:
                            distance = get_distance(a_point, b_point)
                            if distance < min_distance:
                                min_distance = distance
                                best_connection = (a_point, b_point, a_type, b_type, a_line, b_line)
                                best_a_idx = i
                                best_b_idx = j
            
            if best_connection:
                a_point, b_point, a_type, b_type, a_line, b_line = best_connection
                connection_count += 1
                
                # 根据端点类型确定贝塞尔曲线方向
                curve_start = None
                curve_end = None
                connection_valid = False
                
                if a_type == 'a_start' and b_type == 'b_end':
                    # A线起点连B线终点：符合头尾相连原则
                    curve_start = a_point
                    curve_end = b_point
                    connection_valid = True
                    arcpy.AddMessage(u"第{0}次连接：A组线起点连接B组线终点".format(connection_count))
                elif a_type == 'a_end' and b_type == 'b_start':
                    # A线终点连B线起点：符合头尾相连原则
                    curve_start = b_point
                    curve_end = a_point
                    connection_valid = True
                    arcpy.AddMessage(u"第{0}次连接：B组线起点连接A组线终点".format(connection_count))
                else:
                    # 其他情况：强制连接最近的点
                    curve_start = a_point
                    curve_end = b_point
                    connection_valid = True
                    arcpy.AddMessage(u"第{0}次连接：强制连接最近点对".format(connection_count))
                
                if connection_valid and curve_start and curve_end:
                    # 基于折线头部或尾部一小段延长两条线并找到交点作为控制点
                    try:
                        # 根据连接的端点类型判断使用头部还是尾部一小段
                        # 如果连接的是线的起点，使用头部一小段；如果连接的是终点，使用尾部一小段
                        use_end_segment_a = (a_type == 'a_end')  # A线终点连接时使用尾部一小段
                        use_end_segment_b = (b_type == 'b_end')  # B线终点连接时使用尾部一小段
                        
                        extended_a = extend_line_from_last_segment(a_line, length=1000, use_end_segment=use_end_segment_a)
                        extended_b = extend_line_from_last_segment(b_line, length=1000, use_end_segment=use_end_segment_b)
                        
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
                                    arcpy.AddMessage(u"成功创建第{0}条贝塞尔曲线，使用折线最后一小段延长线交点作为控制点".format(connection_count))
                                else:
                                    arcpy.AddWarning(u"第{0}条贝塞尔曲线创建失败".format(connection_count))
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
                                    arcpy.AddMessage(u"第{0}条贝塞尔曲线使用中点作为控制点创建成功".format(connection_count))
                        else:
                            arcpy.AddWarning(u"第{0}次连接：线延长失败".format(connection_count))
                            
                    except Exception as curve_error:
                        arcpy.AddWarning(u"第{0}次连接创建贝塞尔曲线时出错: {1}".format(connection_count, curve_error))
                
                # 移除已连接的线段
                remaining_a.pop(best_a_idx)
                remaining_b.pop(best_b_idx)
                arcpy.AddMessage(u"移除已连接线段，A组剩余{0}条，B组剩余{1}条".format(len(remaining_a), len(remaining_b)))
            else:
                arcpy.AddWarning(u"未找到有效连接，退出循环")
                break
        
        # A组剩下最后一条线时，与B组所有剩余线全连接
        if len(remaining_a) == 1 and len(remaining_b) > 0:
            last_a_line = remaining_a[0]
            arcpy.AddMessage(u"A组剩余最后1条线，开始与B组剩余{0}条线全连接".format(len(remaining_b)))
            
            for i, line_b in enumerate(remaining_b):
                connection_count += 1
                
                # 找到A线和B线的最近点对
                a_points = [
                    (last_a_line.firstPoint, 'a_start'),
                    (last_a_line.lastPoint, 'a_end')
                ]
                b_points = [
                    (line_b.firstPoint, 'b_start'),
                    (line_b.lastPoint, 'b_end')
                ]
                
                # 先找到最近的点对
                min_distance = float('inf')
                closest_a_point = None
                closest_b_point = None
                closest_a_type = None
                closest_b_type = None
                
                for a_point, a_type in a_points:
                    for b_point, b_type in b_points:
                        distance = get_distance(a_point, b_point)
                        if distance < min_distance:
                            min_distance = distance
                            closest_a_point = a_point
                            closest_b_point = b_point
                            closest_a_type = a_type
                            closest_b_type = b_type
                
                # 判断最近点对是否符合头尾相连原则
                best_connection = None
                if (closest_a_type == 'a_start' and closest_b_type == 'b_end') or (closest_a_type == 'a_end' and closest_b_type == 'b_start'):
                    if closest_a_type == 'a_start' and closest_b_type == 'b_end':
                        best_connection = (closest_a_point, closest_b_point, 'a_start_b_end')
                    elif closest_a_type == 'a_end' and closest_b_type == 'b_start':
                        best_connection = (closest_b_point, closest_a_point, 'b_start_a_end')
                
                if best_connection:
                    best_a_point, best_b_point, connection_type = best_connection
                    if connection_type == 'a_start_b_end':
                        arcpy.AddMessage(u"最后阶段第{0}次连接：A组线起点连接B组线{1}终点".format(connection_count, i+1))
                    else:
                        arcpy.AddMessage(u"最后阶段第{0}次连接：A组线终点连接B组线{1}起点".format(connection_count, i+1))
                    try:
                        # 根据连接的端点类型判断使用头部还是尾部一小段
                        use_end_segment_a = (closest_a_type == 'a_end')  # A线终点连接时使用尾部一小段
                        use_end_segment_b = (closest_b_type == 'b_end')  # B线终点连接时使用尾部一小段
                        
                        extended_a = extend_line_from_last_segment(last_a_line, length=1000, use_end_segment=use_end_segment_a)
                        extended_b = extend_line_from_last_segment(line_b, length=1000, use_end_segment=use_end_segment_b)
                        
                        if extended_a and extended_b:
                            intersection = find_line_intersection(extended_a, extended_b)
                            
                            if intersection:
                                curve = create_bezier_curve_with_control_point(
                                    best_a_point, best_b_point, intersection, 
                                    fullness_factor, num_points, spatial_ref
                                )
                            else:
                                # 使用中点作为控制点
                                mid_x = (best_a_point.X + best_b_point.X) / 2
                                mid_y = (best_a_point.Y + best_b_point.Y) / 2
                                mid_point = arcpy.Point(mid_x, mid_y)
                                
                                curve = create_bezier_curve_with_control_point(
                                    best_a_point, best_b_point, mid_point,
                                    fullness_factor, num_points, spatial_ref
                                )
                            
                            if curve:
                                curves.append(curve)
                                arcpy.AddMessage(u"最后阶段：成功创建第{0}条贝塞尔曲线（A线与B组第{1}条线连接）".format(connection_count, i+1))
                            else:
                                arcpy.AddWarning(u"最后阶段：第{0}条贝塞尔曲线创建失败".format(connection_count))
                        else:
                            arcpy.AddWarning(u"最后阶段：线延长失败")
                            
                    except Exception as curve_error:
                        arcpy.AddWarning(u"最后阶段创建贝塞尔曲线时出错: {0}".format(curve_error))
                else:
                    # 如果没有符合头尾相连原则的连接，跳过这条线
                    arcpy.AddWarning(u"最后阶段：A组线与B组第{0}条线没有符合头尾相连原则的连接点，跳过".format(i+1))
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

def create_cross_group_straight_connections(group_c, group_d, spatial_ref):
    """
    在C组和D组线要素之间创建直线连接
    先确定全局连接方向（C头D尾或C尾D头），然后按坐标排序进行连接
    """
    try:
        lines = []
        arcpy.AddMessage(u"开始创建C组({0}条)和D组({1}条)之间的直线连接".format(len(group_c), len(group_d)))
        
        if not group_c or not group_d:
            arcpy.AddWarning(u"C组或D组线要素为空，无法创建连接")
            return lines
        
        # 第一步：确定全局连接方向（C头D尾 或 C尾D头）
        # 随便取第一条C线和第一条D线来判断
        first_c = group_c[0]
        first_d = group_d[0]
        
        c_start = first_c.firstPoint
        c_end = first_c.lastPoint
        d_start = first_d.firstPoint
        d_end = first_d.lastPoint
        
        # 计算两种连接方向的距离
        c_head_d_tail_distance = get_distance(c_start, d_end)  # C头D尾
        c_tail_d_head_distance = get_distance(c_end, d_start)  # C尾D头
        
        # 确定全局连接方向
        if c_head_d_tail_distance <= c_tail_d_head_distance:
            connection_mode = "c_head_d_tail"
            arcpy.AddMessage(u"确定连接方向：C头连D尾（距离：{0:.2f}）".format(c_head_d_tail_distance))
        else:
            connection_mode = "c_tail_d_head"
            arcpy.AddMessage(u"确定连接方向：C尾连D头（距离：{0:.2f}）".format(c_tail_d_head_distance))
        
        # 第二步：根据连接方向对线段进行排序
        # 根据CD中任意一条线的方向，顺时针转90度的方向作为判断点顺序大小
        
        # 选择一条参考线来确定排序方向（优先选择C组第一条线）
        reference_line = group_c[0] if group_c else group_d[0]
        ref_start = reference_line.firstPoint
        ref_end = reference_line.lastPoint
        
        # 计算参考线的方向向量
        line_dx = ref_end.X - ref_start.X
        line_dy = ref_end.Y - ref_start.Y
        
        # 顺时针旋转90度得到垂直方向向量
        # 原向量(dx, dy) 顺时针旋转90度后变为(dy, -dx)
        perp_dx = line_dy
        perp_dy = -line_dx
        
        # 归一化垂直向量（避免除零）
        perp_length = (perp_dx * perp_dx + perp_dy * perp_dy) ** 0.5
        if perp_length > 0:
            perp_dx = perp_dx / perp_length
            perp_dy = perp_dy / perp_length
        else:
            # 如果参考线长度为0，使用默认方向
            perp_dx = 1.0
            perp_dy = 0.0
        
        arcpy.AddMessage(u"参考线方向：({0:.3f}, {1:.3f})，垂直排序方向：({2:.3f}, {3:.3f})".format(
            line_dx, line_dy, perp_dx, perp_dy))
        
        # 根据连接方向确定要排序的端点
        if connection_mode == "c_head_d_tail":
            # C头D尾：分析C的头点和D的尾点坐标
            c_points = [line.firstPoint for line in group_c]
            d_points = [line.lastPoint for line in group_d]
        else:
            # C尾D头：分析C的尾点和D的头点坐标
            c_points = [line.lastPoint for line in group_c]
            d_points = [line.firstPoint for line in group_d]
        
        # 定义排序函数：计算点在垂直方向上的投影
        def get_projection_value(point):
            return point.X * perp_dx + point.Y * perp_dy
        
        # 根据连接方向定义排序键函数
        if connection_mode == "c_head_d_tail":
            def get_c_key(line):
                return get_projection_value(line.firstPoint)
            def get_d_key(line):
                return get_projection_value(line.lastPoint)
        else:
            def get_c_key(line):
                return get_projection_value(line.lastPoint)
            def get_d_key(line):
                return get_projection_value(line.firstPoint)
        
        c_lines_sorted = sorted(group_c, key=get_c_key)
        d_lines_sorted = sorted(group_d, key=get_d_key)
        arcpy.AddMessage(u"按垂直方向投影排序完成")
        
        # 创建剩余线段列表的副本
        remaining_c = list(c_lines_sorted)
        remaining_d = list(d_lines_sorted)
        
        connection_count = 0
        
        # 第三步：按顺序连接，直到C组剩下最后一条线
        while len(remaining_c) > 1 and len(remaining_d) > 0:
            connection_count += 1
            
            # 取C组和D组的第一条线进行连接
            line_c = remaining_c[0]
            line_d = remaining_d[0]
            
            # 根据连接方向确定连接点
            if connection_mode == "c_head_d_tail":
                end_point = line_c.firstPoint
                start_point = line_d.lastPoint
                connection_desc = "D尾连C头"
            else:
                start_point = line_c.lastPoint
                end_point = line_d.firstPoint
                connection_desc = "C尾连D头"
            
            # 创建直线连接
            try:
                line_points = [start_point, end_point]
                straight_line = arcpy.Polyline(arcpy.Array(line_points), spatial_ref)
                lines.append(straight_line)
                distance = get_distance(start_point, end_point)
                arcpy.AddMessage(u"第{0}次连接：{1}（距离：{2:.2f}）".format(connection_count, connection_desc, distance))
            except Exception as line_error:
                arcpy.AddWarning(u"第{0}次连接创建直线时出错: {1}".format(connection_count, line_error))
            
            # 移除已连接的线段
            remaining_c.pop(0)
            remaining_d.pop(0)
            arcpy.AddMessage(u"移除已连接线段，C组剩余{0}条，D组剩余{1}条".format(len(remaining_c), len(remaining_d)))
        
        # 第四步：C组剩下最后一条线时，与D组所有剩余线按顺序连接
        if len(remaining_c) == 1 and len(remaining_d) > 0:
            last_c_line = remaining_c[0]
            arcpy.AddMessage(u"C组剩余最后1条线，开始与D组剩余{0}条线按顺序连接".format(len(remaining_d)))
            
            for i, line_d in enumerate(remaining_d):
                connection_count += 1
                
                # 根据连接方向确定连接点
                if connection_mode == "c_head_d_tail":
                    end_point = last_c_line.firstPoint
                    start_point = line_d.lastPoint
                    connection_desc = "C头连D尾"
                else:
                    start_point = last_c_line.lastPoint
                    end_point = line_d.firstPoint
                    connection_desc = "C尾连D头"
                
                try:
                    # 创建直线连接
                    line_points = [start_point, end_point]
                    straight_line = arcpy.Polyline(arcpy.Array(line_points), spatial_ref)
                    lines.append(straight_line)
                    distance = get_distance(start_point, end_point)
                    arcpy.AddMessage(u"最后阶段第{0}次连接：{1}（距离：{2:.2f}）".format(connection_count, connection_desc, distance))
                except Exception as line_error:
                    arcpy.AddWarning(u"最后阶段创建直线时出错: {0}".format(line_error))
        
        arcpy.AddMessage(u"CD组直线连接完成，成功创建 {0} 条直线".format(len(lines)))
        return lines
        
    except Exception as e:
        arcpy.AddError(u"create_cross_group_straight_connections函数执行失败: {0}".format(e))
        return []

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

        # 设置默认的曲线参数
        fullness_factor = 0.5  # 默认曲线饱满度
        num_points = 20  # 默认曲线平滑度
        arcpy.AddMessage(u"使用默认曲线参数 - 饱满度: {}, 平滑度: {}".format(fullness_factor, num_points))

        arcpy.AddMessage(u"正在获取参数 1 (是否采用曲线连接)...")
        use_curve_connection = arcpy.GetParameter(1)
        if use_curve_connection is None:
            use_curve_connection = True  # 默认勾选，采用曲线连接
            arcpy.AddWarning(u"连接方式参数无效，使用默认值: 曲线连接")
        
        if use_curve_connection:
            connection_mode = u"曲线连接"
        else:
            connection_mode = u"直线连接"
        
        arcpy.AddMessage(u"参数 1 (连接方式): {}".format(connection_mode))

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
            # arcpy.AddMessage(u"【debug1】使用CopyFeatures方法处理选中要素...")
            try:
                # 创建临时要素类来存储选择集，使用系统临时目录
                import tempfile
                temp_dir = tempfile.gettempdir()
                temp_fc = os.path.join(temp_dir, "temp_selected_features.shp")
                # arcpy.AddMessage(u"【debug】创建临时要素类: {0}".format(temp_fc))
                
                # 删除可能存在的临时文件
                if arcpy.Exists(temp_fc):
                    arcpy.Delete_management(temp_fc)
                
                # 使用CopyFeatures复制选择集到临时要素类
                arcpy.CopyFeatures_management(input_features, temp_fc)
                # arcpy.AddMessage(u"【debug】CopyFeatures执行成功")
                
                # 从临时要素类读取要素
                # arcpy.AddMessage(u"【debug】从临时要素类读取要素...")
                with arcpy.da.SearchCursor(temp_fc, ["SHAPE@", "OID@"]) as cursor:
                    # arcpy.AddMessage(u"【debug】SearchCursor创建成功，开始遍历...")
                    for row in cursor:
                        geometry = row[0]
                        oid = row[1]
                        feature_count += 1
                        # # arcpy.AddMessage(u"【debug】读取到要素 OID: {0}".format(oid))
                        
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
                        # arcpy.AddMessage(u"【debug】临时要素类已清理")
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
                # arcpy.AddMessage(u"【debug】fallback使用数据源: {0}".format(data_source))
                # arcpy.AddMessage(u"【debug】尝试使用SearchCursor读取所有要素...")
                with arcpy.da.SearchCursor(data_source, ["SHAPE@", "OID@"]) as cursor:
                    # arcpy.AddMessage(u"【debug】SearchCursor创建成功，开始遍历...")
                    for row in cursor:
                         geometry = row[0]
                         oid = row[1]
                         feature_count += 1
                         # arcpy.AddMessage(u"【debug】读取到要素 OID: {0}".format(oid))
                         
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
                        # arcpy.AddMessage(u"【debug】尝试使用备用数据源: {0}".format(alternative_source))
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

        # 根据连接方式参数选择不同的分组和连接策略
        if use_curve_connection:
            # 曲线连接：使用斜率容差分类线要素为AB两组
            arcpy.AddMessage(u"执行曲线连接方式，开始按斜率分类线要素...")
            
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
        else:
            # 直线连接：使用距离分类线要素为CD两组
            arcpy.AddMessage(u"执行直线连接方式，开始按距离分类线要素...")
            
            try:
                # 按距离分类线要素
                group_c, group_d = classify_lines_by_distance(lines)
                
                if len(group_c) == 0 or len(group_d) == 0:
                    arcpy.AddWarning(u"无法将线要素分为两组，C组: {0} 条，D组: {1} 条".format(len(group_c), len(group_d)))
                    arcpy.AddMessage(u"回退到原有的成对连接方式...")
                    raise Exception(u"线要素分组失败")
                
                arcpy.AddMessage(u"成功分组：C组 {0} 条线，D组 {1} 条线".format(len(group_c), len(group_d)))
                
                # 创建CD组之间的直线连接
                curves = create_cross_group_straight_connections(
                    group_c, group_d, spatial_ref
                )
                
                if not curves or len(curves) == 0:
                    arcpy.AddWarning(u"CD组连接未能创建任何直线，回退到原有连接方式")
                    raise Exception(u"CD组连接创建直线失败")
                    
            except Exception as group_error:
                arcpy.AddWarning(u"CD组连接失败: {0}".format(group_error))
                arcpy.AddMessage(u"回退到原有的成对连接方式...")
        
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
                 connection_type_msg = u"贝塞尔曲线" if use_curve_connection else u"直线连接"
                 arcpy.AddMessage(u"通过临时要素类成功添加了 {0} 条{1}".format(curves_created, connection_type_msg))
                 
                 # 清理临时文件
                 try:
                     if arcpy.Exists(temp_curves_fc):
                         arcpy.Delete_management(temp_curves_fc)
                 except:
                     pass

        final_msg = u"贝塞尔曲线创建成功！" if use_curve_connection else u"直线连接创建成功！"
        arcpy.AddMessage(u"{0}已添加到输入要素类中".format(final_msg))

    except Exception as e:
        arcpy.AddError(u"脚本执行出错: {0}".format(e))
        import traceback
        arcpy.AddError(traceback.format_exc())
