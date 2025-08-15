# -*- coding: utf-8 -*-
import arcpy
import math
import os

def get_distance(p1, p2):
    """计算两个 arcpy.Point 之间的欧氏距离"""
    return math.sqrt((p1.X - p2.X)**2 + (p1.Y - p2.Y)** 2)

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

        # 确定要素ID字段名（解决Shapefile和GDB的差异）
        id_field = "FID" 
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
                # 创建临时要素类来存储选中的要素
                import tempfile
                import os
                temp_dir = tempfile.gettempdir()
                temp_fc = os.path.join(temp_dir, "temp_selected.shp")
                arcpy.management.CopyFeatures(input_features, temp_fc)
                arcpy.AddMessage(u"【debug2】已复制选中要素到临时位置")
                
                # 从临时要素类读取数据
                with arcpy.da.SearchCursor(temp_fc, ["SHAPE@", "OID@"]) as cursor:
                    row_count = 0
                    for row in cursor:
                        row_count += 1
                        geometry = row[0]
                        oid = row[1]
                        feature_count += 1
                        arcpy.AddMessage(u"正在处理选中的要素 OID: {0}".format(oid))
                        
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
                    
                    arcpy.AddMessage(u"【debug3】从临时要素类读取了 {0} 个要素".format(row_count))
                
                # 清理临时要素类
                try:
                    arcpy.management.Delete(temp_fc)
                except:
                    pass  # 忽略删除错误
                    
            except Exception as cursor_error:
                arcpy.AddWarning(u"使用CopyFeatures方法失败: {0}".format(cursor_error))
                arcpy.AddMessage(u" fallback 到处理所有要素...")
                has_user_selection = False  # 失败则处理所有要素

        # 如果没有选择集或选择集处理失败，读取所有要素
        if not has_user_selection or len(lines) == 0:
            arcpy.AddMessage(u"读取所有线要素...")
            try:
                # 使用底层数据源路径
                data_source = desc.catalogPath if hasattr(desc, 'catalogPath') else input_features
                arcpy.AddMessage(u"【debug】fallback使用数据源: {0}".format(data_source))
                with arcpy.da.SearchCursor(data_source, ["SHAPE@", "OID@"]) as cursor:
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
            except Exception as cursor_error:
                arcpy.AddWarning(u"读取要素失败: {0}".format(cursor_error))
        
        arcpy.AddMessage(u"总共处理了 {0} 个要素，其中 {1} 个是有效的线要素。".format(feature_count, len(lines)))

        if not lines:
            if feature_count == 0:
                raise ValueError(u"输入要素集为空。请确保图层包含线要素且有数据。")
            else:
                raise ValueError(u"输入的 {0} 个要素中没有找到有效线要素。".format(feature_count))

        # 检查线要素数量是否为偶数
        if len(lines) % 2 != 0:
            arcpy.AddWarning(u"选择了奇数条线，最后一条线将被忽略。")
            lines = lines[:-1]  # 截断最后一条

        # 创建贝塞尔曲线并插入到输出要素类
        with arcpy.da.InsertCursor(output_fc, ["SHAPE@"]) as insert_cursor:
            curves_created = 0
            # 成对处理线要素
            for i in range(0, len(lines), 2):
                if i + 1 >= len(lines):
                    break  # 确保不越界
                    
                try:
                    line1 = lines[i]
                    line2 = lines[i+1]
                    arcpy.AddMessage(u"正在处理线对 {0}-{1}...".format(i+1, i+2))

                    # 获取线要素端点
                    p1_start = line1.firstPoint
                    p1_end = line1.lastPoint
                    p2_start = line2.firstPoint
                    p2_end = line2.lastPoint
                    
                    # 验证端点有效性
                    points = [p1_start, p1_end, p2_start, p2_end]
                    for j, point in enumerate(points):
                        if not point or not hasattr(point, 'X') or not hasattr(point, 'Y'):
                            raise ValueError(u"线要素 {0} 的端点 {1} 无效".format(i+1 if j < 2 else i+2, j % 2 + 1))
                        
                        if not isinstance(point.X, (int, float)) or not isinstance(point.Y, (int, float)):
                            raise ValueError(u"线要素 {0} 的端点坐标无效: ({1}, {2})".format(i+1 if j < 2 else i+2, point.X, point.Y))

                    # 计算最短连接方式
                    dist_sum1 = get_distance(p1_start, p2_start) + get_distance(p1_end, p2_end)
                    dist_sum2 = get_distance(p1_start, p2_end) + get_distance(p1_end, p2_start)

                    if dist_sum1 < dist_sum2:
                        # 连接方式1: start1->start2, end1->end2
                        curve1 = create_bezier_from_endpoints(p1_start, p2_start, p1_end, p2_end, fullness_factor, num_points, spatial_ref)
                        if curve1:
                            insert_cursor.insertRow([curve1])
                            curves_created += 1
                            arcpy.AddMessage(u"成功创建贝塞尔曲线 {0}".format(curves_created))
                        
                        curve2 = create_bezier_from_endpoints(p1_end, p2_end, p1_start, p2_start, fullness_factor, num_points, spatial_ref)
                        if curve2:
                            insert_cursor.insertRow([curve2])
                            curves_created += 1
                            arcpy.AddMessage(u"成功创建贝塞尔曲线 {0}".format(curves_created))
                    else:
                        # 连接方式2: start1->end2, end1->start2
                        curve1 = create_bezier_from_endpoints(p1_start, p2_end, p1_end, p2_start, fullness_factor, num_points, spatial_ref)
                        if curve1:
                            insert_cursor.insertRow([curve1])
                            curves_created += 1
                            arcpy.AddMessage(u"成功创建贝塞尔曲线 {0}".format(curves_created))
                        
                        curve2 = create_bezier_from_endpoints(p1_end, p2_start, p1_start, p2_end, fullness_factor, num_points, spatial_ref)
                        if curve2:
                            insert_cursor.insertRow([curve2])
                            curves_created += 1
                            arcpy.AddMessage(u"成功创建贝塞尔曲线 {0}".format(curves_created))
                            
                except Exception as pair_error:
                    arcpy.AddWarning(u"处理线对 {0}-{1} 时出错: {2}".format(i+1, i+2, pair_error))
                    continue
            
            arcpy.AddMessage(u"总共成功创建了 {0} 条贝塞尔曲线".format(curves_created))

        arcpy.AddMessage(u"贝塞尔曲线创建成功！输出路径: {0}".format(output_fc))

    except Exception as e:
        arcpy.AddError(u"脚本执行出错: {0}".format(e))
        import traceback
        arcpy.AddError(traceback.format_exc())
