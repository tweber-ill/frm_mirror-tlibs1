#include "../string/string.h"

int main()
{
	std::string str(" \t Test<:>12<3<:>xyz    ");
	std::cout << "|" << str << "|" << std::endl;
	tl::trim(str);
	std::cout << "|" << str << "|" << std::endl;


	std::vector<std::string> vecToks;
	tl::get_tokens<std::string, std::string, std::vector<std::string>>(str, ":", vecToks);
	for(const std::string& strTok : vecToks) std::cout << strTok << std::endl;


	std::vector<std::string> vecToks2;
	tl::get_tokens_seq<std::string, std::string, std::vector<std::string>>(str, "<:>", vecToks2);
	std::cout << std::endl;
	for(const std::string& strTok : vecToks2) std::cout << strTok << std::endl;



	std::string str2("TestSeparatorxyzSEPARATOR123456");

	std::vector<std::string> vecToks3;
	tl::get_tokens_seq<std::string, std::string, std::vector<std::string>>(str2, "Separator", vecToks3, 0);
	std::cout << std::endl;
	for(const std::string& strTok : vecToks3) std::cout << strTok << std::endl;

	return 0;
}
