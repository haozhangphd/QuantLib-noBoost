/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2004, 2005, 2007 StatPro Italia srl

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

/*! \file singleton.hpp
    \brief basic support for the singleton pattern
*/

#ifndef quantlib_singleton_hpp
#define quantlib_singleton_hpp

#include <ql/qldefines.hpp>

#ifdef QL_ENABLE_SINGLETON_THREAD_SAFE_INIT
    #if defined(QL_ENABLE_SESSIONS)
            #warning \
                Thread-safe singleton initialization not supported \
                when sessions are enabled.
    #else
        #include <atomic>
        #include <mutex>
                #warning \
                    Thread-safe singleton initialization \
                    may degrade performances.
        #define QL_SINGLETON_THREAD_SAFE_INIT
    #endif
#endif

#include <ql/types.hpp>
#include <memory>
#include <algorithm>
#if defined(QL_PATCH_MSVC)
    #pragma managed(push, off)
#endif
#if defined(QL_PATCH_MSVC)
    #pragma managed(pop)
#endif
#include <map>


#if (_MANAGED == 1) || (_M_CEE == 1)
// One of the Visual C++ /clr modes. In this case, the global instance
// map must be declared as a static data member of the class.
#define QL_MANAGED 1
#else
// Every other configuration. The map can be declared as a static
// variable inside the creation method.
#define QL_MANAGED 0
#endif

namespace QuantLib {

    #if defined(QL_ENABLE_SESSIONS)
    // definition must be provided by the user
    Integer sessionId();
    #endif

    // this is required on VC++ when CLR support is enabled
    #if defined(QL_PATCH_MSVC)
        #pragma managed(push, off)
    #endif

    //! Basic support for the singleton pattern.
    /*! The typical use of this class is:
        \code
        class Foo : public Singleton<Foo> {
            friend class Singleton<Foo>;
          private:
            Foo() {}
          public:
            ...
        };
        \endcode
        which, albeit sub-optimal, frees one from the concerns of
        creating and managing the unique instance and can serve later
        as a single implemementation point should synchronization
        features be added.

        \ingroup patterns
    */
    template <class T>
    class Singleton {
    #if (QL_MANAGED == 1) && !defined(QL_SINGLETON_THREAD_SAFE_INIT)
      private:
        static std::map<Integer, std::shared_ptr<T> > instances_;
    #endif

    #if defined(QL_SINGLETON_THREAD_SAFE_INIT)
      private:
        static std::atomic<T*> instance_;
        static std::mutex mutex_;
    #endif

      public:
        //! access to the unique instance
        static T& instance();
      protected:
        Singleton() {}
	Singleton( const Singleton& ) = delete;
	Singleton& operator=( const Singleton& ) = delete;
    };

    // static member definitions
    
    #if (QL_MANAGED == 1) && !defined(QL_SINGLETON_THREAD_SAFE_INIT)
      template <class T>
      std::map<Integer, std::shared_ptr<T> > Singleton<T>::instances_;
    #endif

    #if defined(QL_SINGLETON_THREAD_SAFE_INIT) 
    template <class T>  std::atomic<T*> Singleton<T>::instance_;
    template <class T> std::mutex Singleton<T>::mutex_;
    #endif
    
    // template definitions

    template <class T>
    T& Singleton<T>::instance() {

        #if (QL_MANAGED == 0) && !defined(QL_SINGLETON_THREAD_SAFE_INIT)
        static std::map<Integer, std::shared_ptr<T> > instances_;
        #endif

        // thread safe double checked locking pattern with atomic memory calls
        #if defined(QL_SINGLETON_THREAD_SAFE_INIT) 

        T* instance =  instance_.load(std::memory_order_consume);
        
        if (!instance) {
            std::scoped_lock guard(mutex_);
            instance = instance_.load(std::memory_order_consume);
            if (!instance) {
                instance = new T();
                instance_.store(instance, std::memory_order_release);
            }
        }

        #else //this is not thread safe

        #if defined(QL_ENABLE_SESSIONS)
        Integer id = sessionId();
        #else
        Integer id = 0;
        #endif

        std::shared_ptr<T>& instance = instances_[id];
        if (!instance)
            instance = std::shared_ptr<T>(new T);

        #endif

        return *instance;
    }

    // reverts the change above
    #if defined(QL_PATCH_MSVC)
        #pragma managed(pop)
    #endif

}

#undef QL_MANAGED

#endif
