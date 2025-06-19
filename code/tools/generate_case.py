with open("./testcases_lws/scene_with_10000_triangles.txt", "w") as f:
    f.write("""
PerspectiveCamera {
    center 0 2 20
    direction 0 -0.1 -1
    up 0 1 0
    angle 35
    width 400
    height 400
}

Background {
    color 0 0 0 
}

Materials {
    numMaterials 9
    Material { diffuseColor 1 1 1 emission 10 10 10 type_d_rl_rr 1 0 0 }
    Material { diffuseColor 1 0.9 0.8 type_d_rl_rr 1 0 0 }
    Material { diffuseColor 1 0.5 0.5 type_d_rl_rr 1 0 0 }
    Material { diffuseColor 0.5 0.5 1 type_d_rl_rr 1 0 0 }
    Material { diffuseColor 1 0.91 0.81 type_d_rl_rr 1 0 0 }
    Material { diffuseColor 1 0.6 0.6 type_d_rl_rr 1 0 0 shininess 0.5 }
    Material { diffuseColor 1 0.93 0.83 type_d_rl_rr 1 0 0 }
    Material { diffuseColor 0.7 1 0.7 type_d_rl_rr 1 0 0 shininess 0.7 }
    Material { diffuseColor 0.8 0.8 1 type_d_rl_rr 1 0 0 shininess 0.9 }
}

Group {
    numObjects 10093
            
    MaterialIndex 1
    """)

    # 生成 72x72 网格的三角形 (71x71 单元，约 10,082 个三角形)
    step = 5.0 / 71.0  # 10 单位范围，71 个间隔
    for i in range(71):
        for j in range(71):
            x0, z0 = -5 + i * step, -5 + j * step
            x1, z1 = -5 + (i + 1) * step, -5 + j * step
            x2, z2 = -5 + i * step, -5 + (j + 1) * step
            x3, z3 = -5 + (i + 1) * step, -5 + (j + 1) * step
            f.write(f"Triangle {{ vertex0 {x0} 0 {z0} vertex1 {x1} 0 {z1} vertex2 {x2} 0 {z2} }}\n")
            f.write(f"Triangle {{ vertex0 {x1} 0 {z1} vertex1 {x3} 0 {z3} vertex2 {x2} 0 {z2} }}\n")

    # 添加其他对象
    f.write("""
    MaterialIndex 0
    Triangle {
        vertex0 2 3.99 3
        vertex1 2 3.99 -1
        vertex2 -2 3.99 -1
    }
    Triangle {
        vertex0 2 3.99 3
        vertex1 -2 3.99 3
        vertex2 -2 3.99 -1
    }
    MaterialIndex 1
    Plane { normal 0 -1 0 offset -4 }
    MaterialIndex 2
    Plane { normal -1 0 0 offset -5 }
    MaterialIndex 3
    Plane { normal 1 0 0 offset -5 }
    MaterialIndex 4
    Plane { normal 0 0 1 offset -4 }
    Plane { normal 0 0 -1 offset -30 }
    MaterialIndex 5
    Sphere { center -3 -2.5 -1 radius 1.5 }
    MaterialIndex 6
    Plane { normal 0 1 0 offset -4 }
    MaterialIndex 7
    Sphere { center 0 -2.5 1 radius 1.5 }
    MaterialIndex 8
    Sphere { center 3 -2.5 3 radius 1.5 }
}
    """)