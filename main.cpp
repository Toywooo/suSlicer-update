
#include <iostream>
#include <Engines/Utility/suCMDParser.h>
#include "slicer.h"
using namespace std;

void helper()
{
	std::cout << "suSlicer is a mesh slicer developed by RMEC of Shanghai University. \n";
	std::cout << "Usage: \n    suSlicer -i [input stl] -o [output contour]" << std::endl;
	std::cout << "Other parameter include -t -tfile" << std::endl;
	std::cout << std::endl;
	std::cout << "-t Specify thickness of layers[default: 0.2]." << std::endl;
	std::cout << "-r Specify thickness config file" << std::endl;

}

void main(int argc, char* argv[])
{
	suCMDParser option(argv, argc);
	std::string inputFile;
	std::string outputFile;
	float thickness = 0;


	if (option.findOption("-h") || argc == 1)
	{
		helper();
		return;
	}

	////////////////////////////////////////// get parameters��ȡ����
	try {
		//get input file
		if (option.findOption("-i"))//������stl�ļ���û�������ļ����׳��쳣
		{
			inputFile = option.findOptionValue("-i");
			if (option.findOption("-o"))
			{
				outputFile = option.findOptionValue("-o");
			}
			else {
				outputFile = inputFile + ".out";
			}			

		}
		else {
			throw std::exception("No input file!");
		}
		if (option.findOption("-t")) //��ȡ������δ���������׳��쳣
		{
			//In order reduce the complex,  we simplely specify one thickness
			std::string strValue = option.findOptionValue("-t");
			thickness = std::stof(strValue);//stof�������Խ��ַ���ת��������
		}
		else {
			throw std::exception("Please specify a legal thickness!");
		}
	}
	catch (std::exception &e)
	{
		std::cout << e.what() << std::endl;
		//what() ���쳣���ṩ��һ���������������ѱ��������쳣�����ء��⽫�����쳣������ԭ��
		helper();
		return;
	}

	//Pre-processing
	//read stl file                     // 0       1     2
	//compute stl model size(length, weight, height)
	//height = max(z) - min(z)
	//std::vector<float> coords, normals; //(x,y,z,x,y,z,x,y,z....)
	//std::vector<unsigned int> tris, solids;//(0,1,2,2,3,0)������  ������
	Trimesh mesh;
	try {
		mesh.read_file(inputFile.c_str());	
		//c_str��Ϊ����c���Լ��ݣ���c������û��string���ͣ��ʱ���ͨ��string�����ĳ�Ա����c_str()��string ����ת����c�е��ַ�����ʽ��	
	}
	catch (std::exception& e) {
		std::cout << e.what() << std::endl;
		return;
	}

	
	std::vector<float> P;
	
	//size_t n = tris.size() / 3;
	//m.read_triangles(coords, normals, tris);���� ����  ������
	float h = mesh.height();
	float p = 0.3;
	while (p < h)
	{
		P.push_back(p);
		p = p + thickness;
		//P.push_back(p);//�����в��z�����������P��
	}
	//std::vector<std::vector<std::vector<int> > > seg_group = Slicing(n, mesh, P);//
	//for each layer
	//for (int i = 0; i < P.size(); i++)
	//{
	//	std::vector<std::vector<int> > segs;//�߶�
	//	//segs = seg_group[i];
	//	/////////////////////////
	//	std::vector<std::vector<std::vector<float> > > contours;//�պ�ͼ��

	//	int q = segs.size();
		
	Slicing(mesh, P);

		//save contours for each layer.
		//..
		
		//contours = contour_construction(q, segs);
}