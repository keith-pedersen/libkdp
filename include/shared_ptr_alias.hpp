#include <memory>

/*! @file shared_ptr_alias.hpp
 *  @brief This file exposes the std::shared_ptr<U> alias constructor, 
 *  for use by Cython implementation of Vector3, whose T() function 
 *  creates a bound reference to the transverse vector (instead of a copy)
 *  by creating a shared_ptr bound to a Vector2*, but deleted as a Vector3*.
 *  @author Keith Pedersen (Keith.David.Pedersen@gmail.com)
 *  @date 2018
*/
template<typename T, typename U>
std::shared_ptr<T> shared_ptr_alias(std::shared_ptr<U> const& owner, T* pointer)
{
	return std::shared_ptr<T>(owner, pointer);
}
