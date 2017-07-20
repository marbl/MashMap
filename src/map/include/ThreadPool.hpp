/**
 * @file    ThreadPool.hpp
 * @brief   implements master-worker thread pooling for parallelizing read mapping
 *          Borrowed and modified code from Mash (github.com/marbl/Mash)
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef ThreadPool_h
#define ThreadPool_h

#include <pthread.h>

/**
 * @brief     generic thread pooling library using C++ pthreads
 * @details   dispatches input tasks to threads after any of them is available.
 *            maintains a thread safe output queue that guarantees that order of 
 *            output is same as input order
 */
template <class TypeInput, class TypeOutput>
class ThreadPool
{
  private:

    struct OutputQueueNode
    {
        // used to preserve input order when outputting
        
        OutputQueueNode * prev;
        OutputQueueNode * next;
        
        TypeOutput * output;
        bool ready;
    };

    unsigned int threadCount;

    pthread_t* threads;

    std::function<TypeOutput* (TypeInput*)> function;
    TypeInput* inputCurrent;
    OutputQueueNode * outputQueueNodeCurrent;

    pthread_mutex_t* mutexInput;
    pthread_mutex_t* mutexOutput;

    pthread_cond_t* condInput;
    pthread_cond_t* condOutput;

    OutputQueueNode* outputQueueHead;
    OutputQueueNode* outputQueueTail;

    bool finished;

  public:
    
    /* Constructor */
    ThreadPool(std::function<TypeOutput* (TypeInput*)> functionNew, unsigned int threadCountNew)
      :
        threadCount(threadCountNew),
        function(functionNew)
  {
    mutexInput = new pthread_mutex_t();
    mutexOutput = new pthread_mutex_t();

    condInput = new pthread_cond_t();
    condOutput = new pthread_cond_t();

    pthread_mutex_init(mutexInput, NULL);
    pthread_mutex_init(mutexOutput, NULL);

    pthread_cond_init(condInput, NULL);
    pthread_cond_init(condOutput, NULL);

    inputCurrent = 0;

    outputQueueHead = 0;
    outputQueueTail = 0;

    finished = false;

    threads = new pthread_t[threadCount];

    for ( int i = 0; i < threadCount; i++ )
    {
      pthread_create(&threads[i], NULL, &ThreadPool::thread, this);
    }
  }

    /* Destructor */
    ~ThreadPool()
    {
      pthread_mutex_lock(mutexInput);
      finished = true;
      pthread_cond_broadcast(condInput);
      pthread_mutex_unlock(mutexInput);

      for ( int i = 0; i < threadCount; i++ )
      {
        pthread_join(threads[i], NULL);
      }

      delete [] threads;

      while ( outputQueueHead != 0 )
      {
        OutputQueueNode * next = outputQueueHead->next;
        delete outputQueueHead;
        outputQueueHead = next;
      }

      delete mutexInput;
      delete mutexOutput;

      delete condInput;
      delete condOutput;
    }

    /* Check if any thread has placed it's output in the queue */
    bool outputAvailable() const
    {
      bool available;

      pthread_mutex_lock(mutexOutput);
      available = outputQueueHead != 0 && outputQueueHead->ready;
      pthread_mutex_unlock(mutexOutput);

      return available;
    }

    /* Pop the output if available; Calling function is responsible for destructing output object later */
    TypeOutput* popOutputWhenAvailable()
    {
      pthread_mutex_lock(mutexOutput);

      if ( outputQueueHead == 0 )
      {
        std::cerr << "ERROR: waiting for output when no output queued\n";
        pthread_mutex_unlock(mutexOutput);
        return 0;
      }

      while ( ! outputQueueHead->ready )
      {
        pthread_cond_wait(condOutput, mutexOutput);
      }

      TypeOutput * output = outputQueueHead->output;

      OutputQueueNode * next = outputQueueHead->next;

      if ( outputQueueTail == outputQueueHead )
      {
        outputQueueTail = 0;
      }

      delete outputQueueHead;
      outputQueueHead = next;
      pthread_mutex_unlock(mutexOutput);

      return output;
    }

    /* Check if any of the threads is still running */
    bool running() const
    {
      bool running;

      pthread_mutex_lock(mutexOutput);
      running = outputQueueHead != 0;
      pthread_mutex_unlock(mutexOutput);

      return running;
    }

    /* Assign job to next available thread (wait if unavailable), thread will destruct the input */
    void runWhenThreadAvailable(TypeInput * input)
    {
      pthread_mutex_lock(mutexInput);

      while ( inputCurrent != 0 )
      {
        pthread_cond_wait(condInput, mutexInput);
      }

      inputCurrent = input;

      // enqueue output while input locked (to preserve order)
      //
      OutputQueueNode * outputQueueNode = new OutputQueueNode();
      outputQueueNode->next = 0;
      outputQueueNode->ready = false;
      //
      pthread_mutex_lock(mutexOutput);
      //
      if ( outputQueueHead == 0 )
      {
        outputQueueHead = outputQueueNode;
      }
      //
      outputQueueNode->prev = outputQueueTail;
      //
      if ( outputQueueTail != 0 )
      {
        outputQueueTail->next = outputQueueNode;
      }
      //
      outputQueueTail = outputQueueNode;
      //
      pthread_mutex_unlock(mutexOutput);

      outputQueueNodeCurrent = outputQueueNode;

      pthread_mutex_unlock(mutexInput);
      pthread_cond_broadcast(condInput);
    }

  private:

    /* Main function that each thread executes */
    static void * thread(void * arg)
    {
      ThreadPool * threadPool = (ThreadPool *)arg;
      TypeInput * input;
      OutputQueueNode * outputQueueNode;

      while ( ! threadPool->finished )
      {
        // wait for input
        //
        pthread_mutex_lock(threadPool->mutexInput);
        //
        while ( ! threadPool->finished && threadPool->inputCurrent == 0 )
        {
          pthread_cond_wait(threadPool->condInput, threadPool->mutexInput);
        }

        if ( threadPool->finished )
        {
          pthread_mutex_unlock(threadPool->mutexInput);
          return 0;
        }
        //
        input = threadPool->inputCurrent;
        outputQueueNode = threadPool->outputQueueNodeCurrent;
        threadPool->inputCurrent = 0;

        pthread_mutex_unlock(threadPool->mutexInput);

        pthread_cond_broadcast(threadPool->condInput);

        // run function
        //
        outputQueueNode->output = threadPool->function(input);

        delete input;

        // signal output
        //
        outputQueueNode->ready = true;
        //
        pthread_mutex_lock(threadPool->mutexOutput);
        pthread_cond_broadcast(threadPool->condOutput);
        pthread_mutex_unlock(threadPool->mutexOutput);
      }

      return NULL;
    }
};

#endif
