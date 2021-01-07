## suSlicer

suSlicer是一个针对stl文件的切片工具，其重要算法参考[Rodrigo Minetto](https://www.sciencedirect.com/science/article/abs/pii/S0010448517301215)等人的工作。

<p float="left">
<img src="doc/intro.png" alt="introduction"  width="600" />
</p>


### Usage:

suSclicer  [[-i] <stl file>] [[-o] <file path>] [[-t] layer_thickness] <float number>] [[-tfile] < file path>]   

* `-i`:  Input file, only stl format is supported.
* `-o`: Output file,  a contours list, as shown in **Output**. [slc](http://paulbourke.net/dataformats/slc/) format file is a further choice.
* `-t`: The thickness of the first sliced layer.
* `-tfile`: The file that defines thickness of each layers

Example:

```
suSlicer -i small.stl  -t  0.2 
```

Output:
```
Z 0.3 0.3   #(Z z dz)Z is a label, z=z coordinate, dz is the layer thickness
x, y
x,y 
x,y
...
C           #end of a contour
x, y
x,y 
x,y
...
C           #end of a contour
Z 0.5 0.2   #new layer
...
C
...
```

### How to start

#### Compile:

Firstly, we use [CMake](https://cmake.org) to generate visual studio solution, then compiling and linking with visual C++.

#### Coding

代码中包含两个工程：主程序suSlicer,以及单元测试程序 unitTest.
主要功能由定义在suSlicer.h和suSlicer.cpp中的四个函数实现：
1. binary_searching
2. Build_triangle_list
3. ComputeIntersection
4. Slicing
5. contour_construction

#### Data Structure
作为一个学习项目，所有数据结构都由std::vector, std::map两种标准库容器实现。


```c++
class Trimesh
{
public:
	Trimesh() { p_mesh = NULL; }
	void read_file(const char* file_name);
	float height() { return maxP_[2] - minP_[2]; }
	float get_f_minzx(unsigned int fi);
	//Ecapsulate funcs
	const float* V(const size_t vi) { return p_mesh->vrt_coords(vi); }
	const unsigned int* F(const size_t fi) { return p_mesh->tri_corner_inds(fi); }
	const float* N(const size_t ni) { return p_mesh->tri_normal(ni); }
	size_t num_F() { return p_mesh->num_tris(); }

	//IO
	void save_to_ascii(char* filename);
	void save_faces_to_ascii(char* filename, std::vector<unsigned int> faces);
	
public:
	// All members are public for convinient using.
	stl_reader::StlMesh <float, unsigned int> *p_mesh;
	std::vector<float> minP_;
	std::vector<float> maxP_;
};
```
其中stl文件的读写临时采用Sebastian Reiter的[stl_reader](https://github.com/sreiter/stl_reader)，遵循BSD-3-Clause开源协议。stl_reader读入的StlMesh中包含

<p float="left">
<img src="doc/vf.png" alt="verts_face"  width="600" />
</p>

```c++
std::vector<TNumber>	coords;   //x,y,z,x,y,z... 坐标列表
std::vector<TNumber>	normals;  //x,y,z,x,y,z... trangle(face）法向量列表
std::vector<TIndex>		tris;     //t1_id0, t1_id1, t1_id2,t2_id0, t2_id1, t2_id2...三角形三个顶点索引列表
std::vector<TIndex>		solids;   //连续三角面片实体列表
```



