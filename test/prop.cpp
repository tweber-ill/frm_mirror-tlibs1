#include "../file/prop.h"

int main()
{
	tl::Prop<std::string> prop;
	prop.Load("tst.ini");

	prop.Add<std::string>("Sec1/Key3", "Test");

	prop.Save("tst.xml");

	std::cout << prop.Query<std::string>("Sec1/Key1") << std::endl;
	return 0;
}
