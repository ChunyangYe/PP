#include"Parafun.h"
#include<string>
#include<iostream>
#include<fstream>
#include<io.h>
#include<direct.h>
using namespace std;
int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		cerr << "Syntax: " << argv[0] << " <input mesh>" << endl;
		return -1;
	}
	string path_save="result/";
	if (access(path_save.c_str(), 0))
    {
        mkdir(path_save.c_str());
    }
	const string input_mesh = argv[1];

	Parafun parafun_m(input_mesh);

	cout << input_mesh << " parameterization begin ..." << endl;
	parafun_m.run_bpe();
	cout << input_mesh << " parameterization finish !" << endl;

	return 0;
}
