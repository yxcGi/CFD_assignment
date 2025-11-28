#ifndef THREADPOOL_H_
#define THREADPOOL_H_

#include <iostream>
#include <thread>
#include <future>
#include <queue>
#include <atomic>
#include <functional>
#include <unordered_map>
#include <chrono>


static constexpr int DEFAULT_INIT_THREAD_NUM = 4;      // 默认初始线程数量
static constexpr int MAX_TASK_QUEUE_NUM = 100;         // 最大任务队列数量

class Thread
{
    using ThreadFunc = std::function<void()>;
public:
    Thread() = delete;
    Thread(ThreadFunc threadFunc)       // 创建线程的时候传入线程函数
        : id_(initId_++)
        , threadFunc_(threadFunc)
    {}
public:
    // 获取线程id
    int getId()
    {
        return id_;
    }

    // 启动线程
    void start()
    {
        // 执行线程函数
        std::thread t{ threadFunc_ };
        t.detach();
    }

private:
    ThreadFunc threadFunc_;
    int id_;
    static int initId_;
};
inline int Thread::initId_ = 0;


class ThreadPool
{
public:
    ~ThreadPool()
    {
        // 所有线程状态：阻塞 || 等待
        stop_ = true;
        notEmpty_.notify_all();  // 通知所有线程退出
        std::unique_lock<std::mutex> lock(mtx_);
        exit_.wait(lock, [this]() { return taskQue_.empty(); });
        std::cout << "ThreadPool stoped." << std::endl;
    }
public:
    static ThreadPool& getInstance()
    {
        static ThreadPool instance;
        return instance;
    }
public:
    // 提交任务
    template<typename Func, typename... Args>
    auto submitTask(Func&& func, Args&&... args) -> std::future<decltype(func(std::forward<Args>(args)...))>
    {
        if (stop_)  // 若线程池已经停止
        {
            throw std::runtime_error("ThreadPool is stopped.");
        }

        // 加锁操作
        {
            std::unique_lock<std::mutex> lock(mtx_);

            using ReType = decltype(func(std::forward<Args>(args)...));

            // 不满才继续放入, 超时抛出异常
            if (!notFull_.wait_for(
                lock,
                std::chrono::seconds(std::chrono::seconds(1)),
                [this]() { return taskQue_.size() < maxTaskQueNum_; }))
            {
                std::cerr << "the taskQue is full, task will be discarded." << std::endl;
                // 抛出异常
                throw std::runtime_error("Task queue is full.");
            }

            // 打包任务
            auto task =
                std::make_shared<std::packaged_task<ReType()>>(
                    std::bind(std::forward<Func>(func), std::forward<Args>(args)...));

            std::future<ReType> res = task->get_future();     // 获取future返回值

            // 任务通过lambda表达式转化为void()放入队列
            taskQue_.emplace([task]() { (*task)(); });
            notEmpty_.notify_all();
            return res;
        }
    }

    // 线程函数
    void threadFunc()
    {
        while (true)      // 不断从任务队列中拿任务
        {
            std::function<void()> task;
            // 加锁
            {
                std::unique_lock<std::mutex> lock(mtx_);
                // 不空才继续取任务
                notEmpty_.wait(lock, [this]() { return !taskQue_.empty() || stop_; });
                // 情况一：不空 没停 -> 继续取任务
                // 情况二：空了 没停 -> 等待
                // 情况三：不空 停了 -> 继续取任务
                // 情况四：空了 停了 -> 返回
                if (stop_ && taskQue_.empty())  // 如果线程池停止且任务队列为空，退出线程
                {
                    exit_.notify_all();
                    return;
                }

                // 取任务
                task = std::move(taskQue_.front());
                taskQue_.pop();
                idleThreadNum_--;
                notFull_.notify_all();
                if (!taskQue_.empty())      // 拿完如果任务队列仍有剩余，继续取任务
                {
                    notEmpty_.notify_all();
                }
            }       // 拿到任务就可以释放锁了
            task();
            idleThreadNum_++;
        }
    }

    void start(int initThreadNum = DEFAULT_INIT_THREAD_NUM, int maxTaskQueNum = MAX_TASK_QUEUE_NUM)     // 启动线程池并设置线程数量
    {
        initThreadNum_ = initThreadNum;     // 初始化线程数量
        maxTaskQueNum_ = maxTaskQueNum;     // 设置最大任务队列数量
        // 创建线程
        for (int i = 0; i < initThreadNum_; ++i)
        {
            auto thread =
                std::make_unique<Thread>(std::bind(&ThreadPool::threadFunc, this));
            pool_.emplace(thread->getId(), std::move(thread));      // 将线程id和线程放入map中
        }
        // 启动线程
        for (auto& thread : pool_)
        {
            thread.second->start();     // 启动线程
            idleThreadNum_++;
        }
    }

    void stop()
    {
        stop_ = true;
    }

private:
    ThreadPool()
        : idleThreadNum_(0)
        , initThreadNum_(0)
        , maxTaskQueNum_(0)
        , stop_(false)
    {}
    ThreadPool(const ThreadPool&) = delete;
    ThreadPool(ThreadPool&&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;
    ThreadPool& operator=(ThreadPool&&) = delete;


private:
    std::queue<std::function<void()>> taskQue_;                 // 任务队列
    std::unordered_map<int, std::unique_ptr<Thread>> pool_;     // 线程池, 包括记录线程id

    std::mutex mtx_;                    // 互斥锁
    std::condition_variable notFull_;   // 队列不满
    std::condition_variable notEmpty_;  // 队列不空
    std::condition_variable exit_;      // 线程退出

    std::atomic_uint idleThreadNum_;    // 空闲线程数量（可变）
    int initThreadNum_;                 // 初始线程数量
    int maxTaskQueNum_;                 // 最大任务队列数量
    
    std::atomic_bool stop_;      // 线程是否停止
};

#endif // THREADPOOL_H_