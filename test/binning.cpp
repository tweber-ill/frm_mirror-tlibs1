/**
 * reduce binning by 2
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

// clang -DNO_IOSTR -o binning binning.cpp ../log/log.cpp -std=c++11 -lstdc++ -lboost_iostreams -lm

#include <iostream>
#include "../file/loaddat.h"

using t_real = double;


int main(int argc, char** argv)
{
	tl::DatFile<t_real, char> dat;
	dat.SetSeparatorChars({':'});

	if(argc < 4)
	{
		std::cerr << "Usage: <in-file> <out-file> <counter-col> [error-col]" << std::endl;
		return -1;
	}

	const char* pcIn = argv[1];
	const char* pcOut = argv[2];

	int iColCnt = tl::str_to_var<int>(std::string(argv[3]));
	int iColErr = -1;
	if(argc >= 5)
		iColErr = tl::str_to_var<int>(std::string(argv[4]));

	dat.Load(pcIn);
	if(!dat)
	{
		std::cerr << "Cannot open \"" << pcIn << "\" for reading." << std::endl;
		return -1;
	}


	std::cout << "Number of columns: " << dat.GetColumnCount() << std::endl;
	std::cout << "Number of rows: " << dat.GetColumn(0).size() << std::endl;

	const auto& hdr = dat.GetHeader();
	for(const auto& pair : hdr)
		std::cout << pair.first << " = " << pair.second << std::endl;


	std::vector<std::vector<t_real>> vecNew;

	for(std::size_t iCol=0; iCol<dat.GetColumnCount(); ++iCol)
	{
		std::vector<t_real> vecNewCol;
		const auto& vecCol = dat.GetColumn(iCol);

		for(std::size_t iRow=0; iRow+1<vecCol.size(); iRow+=2)
		{
			t_real dDat = t_real();

			if(iCol == iColCnt)
				dDat = vecCol[iRow] + vecCol[iRow+1];
			else if(iCol == iColErr)
				dDat = std::sqrt(vecCol[iRow]*vecCol[iRow] + vecCol[iRow+1]*vecCol[iRow+1]);
			else
				dDat = (vecCol[iRow] + vecCol[iRow+1]) / t_real(2);

			vecNewCol.push_back(dDat);
		}

		vecNew.emplace_back(std::move(vecNewCol));
	}



	std::ofstream ofstrOut(pcOut);
	if(!ofstrOut)
	{
		std::cerr << "Cannot open \"" << pcOut << "\" for writing." << std::endl;
		return -1;
	}

	for(std::size_t iRow=0; iRow < vecNew[0].size(); ++iRow)
	{
		for(std::size_t iCol=0; iCol < vecNew.size(); ++iCol)
			ofstrOut << vecNew[iCol][iRow] << "\t";
		ofstrOut << "\n";
	}

	return 0;
}
