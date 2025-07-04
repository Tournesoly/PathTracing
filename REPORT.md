# 计算机图形学大作业-路径追踪

- 按照大作业推荐实现步骤
- 基于PA1框架


## 基础要求

### 基础要求1 Whitted-Style 光线追踪
<!-- - 材质：折射
- 材质：反射
- Shadow Ray：阴影 -->

光线投射不考虑交点的反射和折射，光线跟踪则考虑，继续跟踪新产生的反射光线和折射光线



### 基础要求2 路径追踪
<!-- - 无穷递归：
    - 支持面光源
    - 终止策略使用 Russian Routelette
- PBR（基于物理的渲染）：
    - 使用蒙特卡洛积分算法计算 Radiance
    - 支持漫反射、理想折射、理想反射材质
    - 支持glossy BRDF（glossy是光泽，BRDF是双向反射分布函数）（如 Cook-Torrance BRDF）
- NEE：对光源采样

前面的光线投射和whitted-style光线追踪表面属性单一，不是反射就是折射，没有漫反射。蒙特卡洛光线跟踪则设定表面属性是混合的，重复发射虚拟光线，叠加混合结果。

- 无穷递归：
    - 支持面光源（whitted-style只支持点光源）
    - 终止策略使用 Russian Routelette（基于概率地终止递归，理论上有概率无穷递归）
- PBR（基于物理的渲染）：
    - 使用蒙特卡洛积分算法计算 Radiance （用来计算颜色）
    - 支持漫反射、理想折射、理想反射材质
    - 支持glossy BRDF（glossy是光泽，BRDF是双向反射分布函数，描述了在点p处，从入射方向进入的光如何被反射到出射方向。它是一个概率密度函数，定义了光线在不同方向上的反射强度。BRDF来决定怎么折射反射漫反射）（如 Cook-Torrance BRDF）
- NEE：对光源采样（每次在点p，随机取一个光源上的点来计算直接光照的贡献，而不是靠随机折射反射漫反射到光源来计算，可以消除直接光照的噪声） -->

总结
- 路径追踪在 Whitted Style 光线追踪的基础上增加了漫反射的间接光照，实现完整的全局光照。
- 它通过BRDF为不同材质（如漫反射、镜面反射、折射、光泽材质）定义光线交互的概率分布
- 使用蒙特卡洛采样根据这些概率追踪光线路径。递归追踪采用 Russian Roulette (RR) 概率终止策略，优化效率并支持理论上的无穷递归。
- 路径追踪还支持面光源，通过随机采样模拟软阴影等真实效果。**NEE（Next Event Estimation）**通过直接对光源采样（发射阴影光线）计算直接光照，显著减少噪声。

### 基础要求3 对比项
- Whitted-Style光线追踪与路径追踪的对比分析
- 无NEE的路径追踪和有NEE的路径追踪的对比分析
- 

## 加分项
TODO
