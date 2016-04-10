// gcc -o eval test/eval.cpp log/log.cpp -std=c++11 -lstdc++ -lm

#include "../string/eval.h"
#include <iostream>

int main()
{
	auto result = tl::eval_expr<std::string, double>("(2^2 - 8) / (-2^1)");
	std::cout << "OK: " << result.first;
	std::cout << ", value: " << result.second << std::endl;

	/*auto result2 = tl::eval_expr<std::wstring, double>(L"2*pi * 100 / 10");
	std::wcout << L"OK: " << result2.first;
	std::wcout << L", value: " << result2.second << std::endl;*/
	return 0;
}
