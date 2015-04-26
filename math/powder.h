/*
 * Powder peaks
 * @author Tobias Weber
 * @date apr-2015
 * @license GPLv2 or GPLv3
 */
 
#ifndef __POWDER_H__
#define __POWDER_H__
 
#include <unordered_set>
#include <tuple>
#include <boost/functional/hash.hpp>
 
namespace tl {

template<class t_int = int>
class Powder
{
	public:
		typedef std::tuple<t_int, t_int, t_int> t_peak;
	
	private:
		static size_t hash_peak(const t_peak& peak)
		{
			boost::hash<t_int> hsh;

			std::size_t iseed = 0;
			boost::hash_combine(iseed, hsh(std::get<0>(peak)));
			boost::hash_combine(iseed, hsh(std::get<1>(peak)));
			boost::hash_combine(iseed, hsh(std::get<2>(peak)));

			return iseed;
		}
		static size_t hash_peak_unique(const t_peak& peak)
		{
			boost::hash<t_int> hsh;
			
			std::size_t iHsh0 = hsh(std::abs(std::get<0>(peak)));
			std::size_t iHsh1 = hsh(std::abs(std::get<1>(peak)));
			std::size_t iHsh2 = hsh(std::abs(std::get<2>(peak)));

			// TODO: find a hash_combine where call order doesn't matter
			std::size_t iseed0 = 0;
			boost::hash_combine(iseed0, iHsh0);
			boost::hash_combine(iseed0, iHsh1);
			boost::hash_combine(iseed0, iHsh2);

			std::size_t iseed1 = 0;
			boost::hash_combine(iseed1, iHsh0);
			boost::hash_combine(iseed1, iHsh2);
			boost::hash_combine(iseed1, iHsh1);

			std::size_t iseed2 = 0;
			boost::hash_combine(iseed2, iHsh1);
			boost::hash_combine(iseed2, iHsh0);
			boost::hash_combine(iseed2, iHsh2);

			std::size_t iseed3 = 0;
			boost::hash_combine(iseed3, iHsh1);
			boost::hash_combine(iseed3, iHsh2);
			boost::hash_combine(iseed3, iHsh0);

			std::size_t iseed4 = 0;
			boost::hash_combine(iseed4, iHsh2);
			boost::hash_combine(iseed4, iHsh0);
			boost::hash_combine(iseed4, iHsh1);

			std::size_t iseed5 = 0;
			boost::hash_combine(iseed5, iHsh2);
			boost::hash_combine(iseed5, iHsh1);
			boost::hash_combine(iseed5, iHsh0);
			
			return iseed0 + iseed1 + iseed2 + iseed3 + iseed4 + iseed5;
		}
		
		static bool equ_peak(const t_peak& peak0, const t_peak& peak1)
		{
			return std::get<0>(peak0) == std::get<0>(peak1) &&
				std::get<1>(peak0) == std::get<1>(peak1) &&
				std::get<2>(peak0) == std::get<2>(peak1);
		}
		static bool equ_peak_unique(const t_peak& peak0, const t_peak& peak1)
		{
			return hash_peak_unique(peak0)==hash_peak_unique(peak1);
		}
		
	public:
		typedef std::unordered_set<t_peak, decltype(&hash_peak), decltype(&equ_peak)> t_peaks;
		typedef std::unordered_set<t_peak, decltype(&hash_peak_unique), decltype(&equ_peak_unique)> t_peaks_unique;

	protected:
		t_peaks m_peaks;
		t_peaks_unique m_peaks_unique;
	
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
		
		void clear()
		{
			m_peaks.clear();
			m_peaks_unique.clear();
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
