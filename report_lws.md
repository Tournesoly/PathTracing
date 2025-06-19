- main 通过argv[3] 控制tracer
- 在tracer.hpp 中实现三种tracer

- 首先 whitted
    - 漫反射
        - 阴影（注意传入引用，咔死我了
    - 反射 注意抬离
    - 折射 注意抬离

- 然后 simplepathtracing
    - 无穷递归
        - 漫反射
        - 反射
        - 折射
        - 俄罗斯轮盘赌
    - 面光源（自发光材质
- neepathracing
- glossypathpacing
- bvhpathpacing

