#include <memory>

template<typename T, typename U>
std::shared_ptr<T> shared_ptr_alias(std::shared_ptr<U> const& owner, T* pointer)
{
	return std::shared_ptr<T>(owner, pointer);
}
