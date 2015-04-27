/*
 * Powder peaks
 * @author Tobias Weber
 * @date apr-2015
 * @license GPLv2 or GPLv3
 */

#ifndef __POWDER_H__
#define __POWDER_H__

#include <tuple>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <initializer_list>
#include "lattice.h"
#include "../helper/hash.h"

namespace tl {

template<class t_int=int, class t_real=double>
class Powder
{
	public:
		typedef std::tuple<t_int, t_int, t_int> t_peak;

	private:
		static size_t hash_peak(const t_peak& peak)
		{
			return tl::hash_ordered<std::initializer_list<t_int>>
				({
					std::get<0>(peak),
					std::get<1>(peak),
					std::get<2>(peak)
				});
		}
		static size_t hash_peak_unique(const t_peak& peak)
		{
			return tl::hash_unordered<std::initializer_list<t_int>>
				({
					std::abs(std::get<0>(peak)),
					std::abs(std::get<1>(peak)),
					std::abs(std::get<2>(peak))
				});
		}

		static bool equ_peak(const t_peak& peak0, const t_peak& peak1)
		{
			return std::get<0>(peak0) == std::get<0>(peak1) &&
				std::get<1>(peak0) == std::get<1>(peak1) &&
				std::get<2>(peak0) == std::get<2>(peak1);
		}
		static bool equ_peak_unique(const t_peak& peak0, const t_peak& peak1)
		{
			std::vector<t_int> pk0({ std::get<0>(peak0), std::get<1>(peak0), std::get<2>(peak0) });
			std::vector<t_int> pk1({ std::get<0>(peak1), std::get<1>(peak1), std::get<2>(peak1) });

			std::transform(pk0.begin(), pk0.end(), pk0.begin(), (t_int(*)(t_int))std::abs);
			std::transform(pk1.begin(), pk1.end(), pk1.begin(), (t_int(*)(t_int))std::abs);

			std::sort(pk0.begin(), pk0.end());
			std::sort(pk1.begin(), pk1.end());

			return pk0 == pk1;
		}

	public:
		typedef std::unordered_set<t_peak, decltype(&hash_peak), decltype(&equ_peak)> t_peaks;
		typedef std::unordered_set<t_peak, decltype(&hash_peak_unique), decltype(&equ_peak_unique)> t_peaks_unique;

	protected:
		t_peaks m_peaks;
		t_peaks_unique m_peaks_unique;

		// associated reciprocal lattice
		const Lattice<t_real> *m_pLatticeRecip = nullptr;

	public:
		Powder()
			: m_peaks(10, &hash_peak, &equ_peak),
			  m_peaks_unique(10, &hash_peak_unique, &equ_peak_unique)
		{}

		void AddPeak(t_int h, t_int k, t_int l)
		{
			t_peak peak(h,k,l);
			t_peak peakabs(std::abs(h),std::abs(k),std::abs(l));

			m_peaks.insert(peak);
			m_peaks_unique.insert(peakabs);
		}

		const t_peaks& GetPeaks() const { return m_peaks; }
		const t_peaks& GetUniquePeaks() const { return m_peaks_unique; }

		bool HasPeak(int h, int k, int l) const
		{
			t_peak peak(h,k,l);
			return m_peaks.find(peak) != m_peaks.end();
		}
		bool HasUniquePeak(int h, int k, int l) const
		{
			t_peak peak(h,k,l);
			return m_peaks_unique.find(peak) != m_peaks_unique.end();
		}

		std::size_t GetMultiplicity(t_int h, t_int k, t_int l) const
		{
			t_peak peak(h,k,l);

			std::size_t iMult = 0;
			for(const t_peak& pk : m_peaks)
			{
				if(equ_peak_unique(pk, peak))
					++iMult;
			}

			return iMult;
		}

		void SetRecipLattice(const Lattice<t_real>* pLatt)
		{
			m_pLatticeRecip = pLatt;
		}

		ublas::vector<t_real> GetRecipLatticePos(double dh, double dk, double dl) const
		{
			if(m_pLatticeRecip)
				return m_pLatticeRecip->GetPos(dh, dk, dl);

			return ublas::vector<t_real>();
		}

		void clear()
		{
			m_peaks.clear();
			m_peaks_unique.clear();
			m_pLatticeRecip = nullptr;
		}
};

}

#include <ostream>

template<class t_int=int>
std::ostream& operator<<(std::ostream& ostr, const tl::Powder<t_int>& powder)
{
	const typename tl::Powder<t_int>::t_peaks& peaks = powder.GetPeaks();
	const typename tl::Powder<t_int>::t_peaks_unique& peaks_unique = powder.GetUniquePeaks();

	ostr << "Peaks:\n";
	for(const typename tl::Powder<t_int>::t_peak& pk : peaks)
	{
		t_int h = std::get<0>(pk);
		t_int k = std::get<1>(pk);
		t_int l = std::get<2>(pk);

		ostr << "\t(" << h << k << l << ")\n";
	}

	ostr << "Unique Peaks:\n";
	for(const typename tl::Powder<t_int>::t_peak& pk : peaks_unique)
	{
		t_int h = std::get<0>(pk);
		t_int k = std::get<1>(pk);
		t_int l = std::get<2>(pk);

		ostr << "\t(" << h << k << l << ")";
		ostr << ", multiplicity: " << powder.GetMultiplicity(h,k,l) << "\n";
	}

	return ostr;
}

#endif
