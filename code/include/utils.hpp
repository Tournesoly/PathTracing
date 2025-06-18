#ifndef UTILS_HPP
#define UTILS_HPP 

#include <cstdlib>
#include <ctime>


// 初始化随机种子（只需在程序启动时调用一次）
static bool seedInitialized = false;
inline void initRand() {
    if (!seedInitialized) {
        srand(static_cast<unsigned int>(time(0))); // 使用当前时间作为种子
        seedInitialized = true;
    }
}

// 随机数生成函数
inline double RandNum() {
    initRand(); // 确保种子初始化
    return static_cast<double>(rand()) / RAND_MAX; // 返回 [0, 1) 的随机数
}



#endif