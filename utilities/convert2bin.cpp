# include <iostream>
# include <fstream>
# include <sstream>
# include <vector>
# include <map>
# include <cstdlib>
# include <cmath>

using namespace std;

int main(int argc, char **argv)
{
	string line;
	ofstream fpo(argv[1]);
	while(getline(cin, line))
	{
		istringstream cl(line);
		float num;
		while(cl >> num)
			fpo.write(reinterpret_cast<const char*>(&num), sizeof num);
	}
}
