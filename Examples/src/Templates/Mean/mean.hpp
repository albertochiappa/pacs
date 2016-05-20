#ifndef _MEANEXAMPLE_H
#define _MEANEXAMPLE_H
#include <type_traits>
namespace Utility{
  //! A simple template for the mean
  /*!
    I use constexpr so it is more efficient if I pass literals
  */
  template <class T>
  constexpr T  mean (T const & a, T const & b)
  {
    return  0.5*(a+b);
  }
  //! Specialization for pointers
  template <class T>
  constexpr T mean (T* const & a, T* const & b)
  {
    return  0.5*(*a+*b);
  }
  
  //! Example use of variadic templates (only for nerds!)
  /*! 
    A function that computes the sum. First the specialization for a 
    single argument
  */
  template <typename T>
  T Sum(const T & a)
  {
    return a;
  }
  //! Full template for the sum
  template <typename T, typename... Ts>
  typename std::common_type<T,Ts...>::type Sum(const T& a, const Ts&... Args)
  {
    return a + Sum(Args...);
  }

}
#endif
