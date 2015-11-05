/*
 * array helpers
 *
 * @author: tweber
 * @date: nov-2015
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_ARRAYS_H__
#define __TLIBS_ARRAYS_H__

namespace tl {


/**
 * Minimalistic wrapper for plain old arrays
 */
template<class T>
class wrapper_array
{
	public:
		using value_type = T;
		using size_type = std::size_t;

		using iterator = T*;
		using const_iterator = const T*;

		using reference = T&;
		using const_reference = const T&;


	protected:
		T* m_pt = nullptr;
		size_type m_len = 0;

	public:
		wrapper_array(T* p, size_type len)
			: m_pt(p), m_len(len)
		{}

		size_type size() const { return m_len; }

		iterator begin() { return m_pt; }
		iterator end() { return m_pt+m_len; }
		const_iterator begin() const { return m_pt; }
		const_iterator end() const { return m_pt+m_len; }

		T& operator[](size_type i) { return m_pt[i]; }
		const T& operator[](size_type i) const { return m_pt[i]; }
};


}
#endif
