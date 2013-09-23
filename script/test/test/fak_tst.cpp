#include <iostream>

int fak(int a)
{
	if(a == 0 || a == 1)
		return 1;
	else
		return a*fak(a-1);
}

int fib(int a)
{
	if(a == 0 || a == 1)
		return 1;

	return fib(a-1) + fib(a-2);
}


int main()
{
	int a = 0;
	while(a < 10)
	{
		int result = fak(a);
		std::cout << a << "! = " << result << std::endl;

		a = a+1;
	}


	std::cout << "\n\n";

	a = 0;
	while(a < 30)
	{
		int result = fib(a);
		std::cout << "fib(" << a << ") = " << result << std::endl;
		++a;
	}

	return 0;
}
