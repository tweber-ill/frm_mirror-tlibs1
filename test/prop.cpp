#include "../file/prop.h"

int main()
{
	tl::Prop<std::string> prop;
	//prop.Load("tst.ini");
	prop.Add<std::string>("Sec1/Key1", "123");
	prop.Add<std::string>("Sec1/Key2", "Test");
	
	std::map<std::string, std::string> map;
	map["Sec2/Key1"] = "abcde";
	map["Sec2/Key2"] = "xyz";
	prop.Add(map);

	prop.Save("tst.xml");

	std::cout << prop.Query<std::string>("Sec1/Key1") << std::endl;
	return 0;
}
