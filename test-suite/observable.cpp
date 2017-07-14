/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015 Klaus Spanderen

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#include "utilities.hpp"
#include <ql/patterns/observable.hpp>
#include <ql/quotes/simplequote.hpp>

using namespace QuantLib;


namespace {
    class UpdateCounter : public Observer {
      public:
        UpdateCounter() : counter_(0) {}
        void update() {
            ++counter_;
        }
        Size counter() { return counter_; }
      private:
        Size counter_;
    };
}

TEST_CASE("Observable_ObservableSettings", "[Observable]") {

    INFO("Testing observable settings...");

    const std::shared_ptr<SimpleQuote> quote = std::make_shared<SimpleQuote>(100.0);
    UpdateCounter updateCounter;

    updateCounter.registerWith(quote);
    if (updateCounter.counter() != 0) {
        FAIL("update counter value is not zero");
    }

   quote->setValue(1.0);
   if (updateCounter.counter() != 1) {
       FAIL("update counter value is not one");
   }

   ObservableSettings::instance().disableUpdates(false);
   quote->setValue(2.0);
   if (updateCounter.counter() != 1) {
       FAIL("update counter value is not one");
   }
   ObservableSettings::instance().enableUpdates();
   if (updateCounter.counter() != 1) {
       FAIL("update counter value is not one");
   }

   ObservableSettings::instance().disableUpdates(true);
   quote->setValue(3.0);
   if (updateCounter.counter() != 1) {
       FAIL("update counter value is not one");
   }
   ObservableSettings::instance().enableUpdates();
   if (updateCounter.counter() != 2) {
       FAIL("update counter value is not two");
   }

   UpdateCounter updateCounter2;
   updateCounter2.registerWith(quote);
   ObservableSettings::instance().disableUpdates(true);
   for (Size i=0; i < 10; ++i) {
       quote->setValue(Real(i));
   }
   if (updateCounter.counter() != 2) {
       FAIL("update counter value is not two");
   }
   ObservableSettings::instance().enableUpdates();
   if (updateCounter.counter() != 3 || updateCounter2.counter() != 1) {
       FAIL("update counter values are not correct");
   }
}


#ifdef QL_ENABLE_THREAD_SAFE_OBSERVER_PATTERN

#include <list>
#include <thread>

namespace {

    class MTUpdateCounter : public Observer {
      public:
        MTUpdateCounter() : counter_(0) {
            ++instanceCounter_;
        }
        ~MTUpdateCounter() {
            --instanceCounter_;
        }
        void update() {
            ++counter_;
        }
        int counter() { return counter_; }
        static int instanceCounter() { return instanceCounter_; }

      private:
        std::atomic<int> counter_;
        static std::atomic<int> instanceCounter_;
    };

    std::atomic<int> MTUpdateCounter::instanceCounter_(0);

    class GarbageCollector {
      public:
        GarbageCollector() : terminate_(false) { }

        void addObj(const std::shared_ptr<MTUpdateCounter>& updateCounter) {
            std::scoped_lock<std::mutex> lock(mutex_);
            objList.emplace_back(updateCounter);
        }

        void run() {
            while(!terminate_) {
                Size objListSize;
                {
                    std::scoped_lock<std::mutex> lock(mutex_);
                    objListSize = objList.size();
                }

                if (objListSize > 20) {
                    // trigger gc
                    while (objListSize > 0) {
                        std::scoped_lock<std::mutex> lock(mutex_);
                        objList.pop_front();
                        objListSize = objList.size();
                    }
                }

                std::this_thread::sleep_for(std::chrono::milliseconds(2));
            }
            objList.clear();
        }

        void terminate() {
            terminate_ = true;
        }
      private:
        std::mutex mutex_;
        std::atomic<bool> terminate_;

        std::list<std::shared_ptr<MTUpdateCounter> > objList;
    };
}

TEST_CASE("Observable_AsyncGarbagCollector", "[.]") {

    INFO("Testing observer pattern with an asynchronous "
                       "garbage collector (JVM/.NET use case)...");

    // This test core dumps if used with the ordinary implementation
    // of the observer pattern (comparable situation
    // in JVM or .NET eco systems).

    const std::shared_ptr<SimpleQuote> quote = std::make_shared<SimpleQuote>(-1.0);

    GarbageCollector gc;
    std::thread workerThread(&GarbageCollector::run, &gc);

    for (Size i=0; i < 10000; ++i) {
        const std::shared_ptr<MTUpdateCounter> observer = std::make_shared<MTUpdateCounter>();
        observer->registerWith(quote);
        gc.addObj(observer);

        for (Size j=0; j < 10; ++j)
            quote->setValue(Real(j));
    }

    gc.terminate();
    workerThread.join();

    if (MTUpdateCounter::instanceCounter() != 0) {
        FAIL("garbage collection does not work.");
    }
}


TEST_CASE("Observable_MultiThreadingGlobalSettings", "[.]") {
	INFO("Testing observer global settings in a "
		               "multithreading environment...");
	
	const std::shared_ptr<SimpleQuote> quote = std::make_shared<SimpleQuote>(-1.0);

    ObservableSettings::instance().disableUpdates(true);

    GarbageCollector gc;
    std::thread workerThread(&GarbageCollector::run, &gc);

    typedef std::list<std::shared_ptr<MTUpdateCounter> > local_list_type;
    local_list_type localList;

    for (Size i=0; i < 4000; ++i) {
        const std::shared_ptr<MTUpdateCounter> observer = std::make_shared<MTUpdateCounter>();
        observer->registerWith(quote);

        if ((i%4) == 0) {
            localList.emplace_back(observer);
            for (Size j=0; j < 5; ++j)
                quote->setValue(Real(j));
        }
        gc.addObj(observer);
    }

    gc.terminate();
    workerThread.join();

    if (localList.size() != MTUpdateCounter::instanceCounter()) {
        FAIL("garbage collection does not work.");
    }

    for (local_list_type::iterator iter = localList.begin();
        iter != localList.end(); ++iter) {
        if ((*iter)->counter() != 0) {
            FAIL("notification should have been blocked");
        }
    }

    ObservableSettings::instance().enableUpdates();

    for (local_list_type::iterator iter = localList.begin();
        iter != localList.end(); ++iter) {
        if ((*iter)->counter() != 1) {
            FAIL("only one notification should have been sent");
        }
    }
}
#endif
