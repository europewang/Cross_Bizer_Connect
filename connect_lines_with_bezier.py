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
    """根据斜率容差将线段分为两组"""
    if len(lines) < 2:
        return lines, []
    
    tolerance_radians = math.radians(tolerance_degrees)
    angles = [get_line_angle(line) for line in lines]
    
    # 使用第一条线作为参考
    reference_angle = angles[0]
    group_a = [lines[0]]
    group_b = []
    
    for i in range(1, len(lines)):
        angle_diff = angle_difference(angles[i], reference_angle)
        if angle_diff <= tolerance_radians:
            group_a.append(lines[i])
        else:
            group_b.append(lines[i])
    
    # 如果group_b为空，尝试用不同的分组策略
    if not group_b and len(lines) > 1:
        # 找到与第一条线角度差最大的线作为第二组的参考
        max_diff = 0
        max_index = 1
        for i in range(1, len(lines)):
            diff = angle_difference(angles[i], reference_angle)
            if diff > max_diff:
                max_diff = diff
                max_index = i
        
        # 重新分组
        second_reference = angles[max_index]
        group_a = []
        group_b = []
        
        for i, line in enumerate(lines):
            diff_to_first = angle_difference(angles[i], reference_angle)
            diff_to_second = angle_difference(angles[i], second_reference)
            
            if diff_to_first <= diff_to_second:
                group_a.append(line)
            else:
                group_b.append(line)
    
    return group_a, group_b

def extend_line(line, length=1000):
    """延长线段"""
    start_point = line.firstPoint
    end_point = line.lastPoint
    
    # 计算方向向量
    dx = end_point.X - start_point.X
    dy = end_point.Y - start_point.Y
    line_length = math.sqrt(dx*dx + dy*dy)
    
    if line_length == 0:
        return None
    
    # 单位方向向量
    unit_dx = dx / line_length
    unit_dy = dy / line_length
    
    # 向两端延长
    new_start_x = start_point.X - unit_dx * length
    new_start_y = start_point.Y - unit_dy * length
    new_end_x = end_point.X + unit_dx * length
    new_end_y = end_point.Y + unit_dy * length
    
    return ((new_start_x, new_start_y), (new_end_x, new_end_y))

def find_line_intersection(line1_extended, line2_extended):
    """计算两条延长线的交点"""
    x1, y1 = line1_extended[0]
    x2, y2 = line1_extended[1]
    x3, y3 = line2_extended[0]
    x4, y4 = line2_extended[1]
    
    # 计算交点
    denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    if abs(denom) < 1e-10:  # 平行线
        return None
    
    t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom
    
    intersection_x = x1 + t * (x2 - x1)
    intersection_y = y1 + t * (y2 - y1)
    
    return arcpy.Point(intersection_x, intersection_y)

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
            x = ((1-t)**2 * start_point.X + 
                 2*(1-t)*t * control_point.X + 
                 t**2 * end_point.X)
            
            y = ((1-t)**2 * start_point.Y + 
                 2*(1-t)*t * control_point.Y + 
                 t**2 * end_point.Y)
            
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

        arcpy.AddMessage(u"正在获取参数 1 (输出路径)...")
        output_fc = arcpy.GetParameterAsText(1)
        if not output_fc:
            raise ValueError(u"输出要素类路径 (参数 1) 未提供。")
        arcpy.AddMessage(u"参数 1 (输出路径): {0}".format(output_fc))

        arcpy.AddMessage(u"正在获取参数 2 (曲线饱满度)...")
        fullness_factor = arcpy.GetParameter(2)
        if fullness_factor is None:
            raise ValueError(u"曲线饱满度 (参数 2) 未提供。")
        arcpy.AddMessage(u"参数 2 (曲线饱满度): {}".format(fullness_factor))

        arcpy.AddMessage(u"正在获取参数 3 (曲线平滑度)...")
        num_points = arcpy.GetParameter(3)
        if num_points is None or num_points < 2:
            num_points = 20  # 默认值
            arcpy.AddWarning(u"曲线平滑度参数无效，使用默认值: {}".format(num_points))
        arcpy.AddMessage(u"参数 3 (曲线平滑度): {}".format(num_points))

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

        # 创建输出要素类
        arcpy.AddMessage(u"正在创建输出要素类...")
        arcpy.CreateFeatureclass_management(
            out_path=os.path.dirname(output_fc),
            out_name=os.path.basename(output_fc),
            geometry_type="POLYLINE",
            spatial_reference=spatial_ref
        )
        arcpy.AddMessage(u"输出要素类创建成功。")

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

        # 按斜率分组线段
        arcpy.AddMessage(u"正在按斜率分组线段...")
        group_a, group_b = classify_lines_by_slope(lines, tolerance_degrees=30)
        
        arcpy.AddMessage(u"分组结果: A组 {0} 条线，B组 {1} 条线".format(len(group_a), len(group_b)))
        
        if not group_a or not group_b:
            arcpy.AddWarning(u"无法将线段分为两组，可能所有线段都是平行的。")
            arcpy.AddWarning(u"将使用原有的成对连接方式...")
            
            # 回退到原有逻辑
            if len(lines) % 2 != 0:
                arcpy.AddWarning(u"选择了奇数条线，最后一条线将被忽略。")
                lines = lines[:-1]
            
            with arcpy.da.InsertCursor(output_fc, ["SHAPE@"]) as insert_cursor:
                curves_created = 0
                for i in range(0, len(lines), 2):
                    if i + 1 >= len(lines):
                        break
                    try:
                        line1 = lines[i]
                        line2 = lines[i+1]
                        p1_start = line1.firstPoint
                        p1_end = line1.lastPoint
                        p2_start = line2.firstPoint
                        p2_end = line2.lastPoint
                        
                        dist_sum1 = get_distance(p1_start, p2_start) + get_distance(p1_end, p2_end)
                        dist_sum2 = get_distance(p1_start, p2_end) + get_distance(p1_end, p2_start)
                        
                        if dist_sum1 < dist_sum2:
                            curve1 = create_bezier_from_endpoints(p1_start, p2_start, p1_end, p2_end, fullness_factor, num_points, spatial_ref)
                            if curve1:
                                insert_cursor.insertRow([curve1])
                                curves_created += 1
                            curve2 = create_bezier_from_endpoints(p1_end, p2_end, p1_start, p2_start, fullness_factor, num_points, spatial_ref)
                            if curve2:
                                insert_cursor.insertRow([curve2])
                                curves_created += 1
                        else:
                            curve1 = create_bezier_from_endpoints(p1_start, p2_end, p1_end, p2_start, fullness_factor, num_points, spatial_ref)
                            if curve1:
                                insert_cursor.insertRow([curve1])
                                curves_created += 1
                            curve2 = create_bezier_from_endpoints(p1_end, p2_start, p1_start, p2_end, fullness_factor, num_points, spatial_ref)
                            if curve2:
                                insert_cursor.insertRow([curve2])
                                curves_created += 1
                    except Exception as pair_error:
                        arcpy.AddWarning(u"处理线对时出错: {0}".format(pair_error))
                        continue
                arcpy.AddMessage(u"总共成功创建了 {0} 条贝塞尔曲线".format(curves_created))
        else:
            # 新的交叉连接逻辑
            with arcpy.da.InsertCursor(output_fc, ["SHAPE@"]) as insert_cursor:
                curves_created = 0
                
                # 为每条A组线段与每条B组线段创建连接
                for i, line_a in enumerate(group_a):
                    for j, line_b in enumerate(group_b):
                        try:
                            arcpy.AddMessage(u"正在处理A组线段{0}与B组线段{1}的连接...".format(i+1, j+1))
                            
                            # 延长两条线段
                            extended_a = extend_line(line_a)
                            extended_b = extend_line(line_b)
                            
                            if not extended_a or not extended_b:
                                arcpy.AddWarning(u"无法延长线段A{0}或B{1}".format(i+1, j+1))
                                continue
                            
                            # 计算延长线交点作为控制点
                            intersection = find_line_intersection(extended_a, extended_b)
                            
                            if not intersection:
                                arcpy.AddWarning(u"线段A{0}与B{1}的延长线平行，无法计算交点".format(i+1, j+1))
                                continue
                            
                            # 获取线段端点
                            a_start = line_a.firstPoint
                            a_end = line_a.lastPoint
                            b_start = line_b.firstPoint
                            b_end = line_b.lastPoint
                            
                            # 计算四种可能的连接方式，选择距离最短的
                            connections = [
                                (a_start, b_start),
                                (a_start, b_end),
                                (a_end, b_start),
                                (a_end, b_end)
                            ]
                            
                            min_distance = float('inf')
                            best_connection = None
                            
                            for start_point, end_point in connections:
                                distance = get_distance(start_point, end_point)
                                if distance < min_distance:
                                    min_distance = distance
                                    best_connection = (start_point, end_point)
                            
                            if best_connection:
                                start_point, end_point = best_connection
                                
                                # 创建贝塞尔曲线，使用交点作为控制点的参考
                                curve = create_bezier_curve_with_control_point(
                                    start_point, end_point, intersection, 
                                    fullness_factor, num_points, spatial_ref
                                )
                                
                                if curve:
                                    insert_cursor.insertRow([curve])
                                    curves_created += 1
                                    arcpy.AddMessage(u"成功创建A{0}-B{1}连接曲线".format(i+1, j+1))
                                
                        except Exception as connection_error:
                            arcpy.AddWarning(u"处理A{0}-B{1}连接时出错: {2}".format(i+1, j+1, connection_error))
                            continue
                
                arcpy.AddMessage(u"总共成功创建了 {0} 条交叉连接贝塞尔曲线".format(curves_created))

        arcpy.AddMessage(u"贝塞尔曲线创建成功！输出路径: {0}".format(output_fc))

    except Exception as e:
        arcpy.AddError(u"脚本执行出错: {0}".format(e))
        import traceback
        arcpy.AddError(traceback.format_exc())
