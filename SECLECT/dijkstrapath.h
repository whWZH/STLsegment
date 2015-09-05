#ifndef _DIJKSTRA_PATH_FCLLIB_H_
#define _DIJKSTRA_PATH_FCLLIB_H_
#include <queue>
#include <map>
#include <stack>
#include <vector>
#include "..\trimesh\include\TriMesh.h"
#include <utility>
using std::priority_queue;
using std::map;
using std::vector;
using std::stack;

#define MAXNUM 100000.0
//this algorithm instance was only written for TriMesh class .. 
class FuDijkstraNode
{
public:
	FuDijkstraNode():id(-1),pre(-1),dist(MAXNUM){}
	FuDijkstraNode(int i,int p,float d):id(i),pre(p),dist(d){}
	FuDijkstraNode(const FuDijkstraNode& rhs)
	{
		id   = rhs.id;
		pre  = rhs.pre;
		dist = rhs.dist;
	}

	FuDijkstraNode& operator= (const FuDijkstraNode& rhs)
	{
		if (this != &rhs)
		{
			id   = rhs.id;
			pre  = rhs.pre;
			dist = rhs.dist;
		}
		return *this;
	}

	//���ȶ����ǰ�����Ԫ�ط���ǰ�棬��������Ҫ����С�ķ���
    //ǰ�棬�������� '< 'operator ������if a > b ,����
    //true��������' < 'operator ���������Ԫ�طź���
	friend bool operator< (const FuDijkstraNode& n1,const FuDijkstraNode& n2)
	{
		return n1.dist > n2.dist;     
	}
                                  
	int    id;
	int    pre;
	float dist;
};

class DijkstraPath
{
public:
	DijkstraPath():themesh(NULL),StartNode(-1),EndNode(-1){}
	DijkstraPath(TriMesh* mesh,int S,int E):themesh(mesh),
		StartNode(S),EndNode(E){}
	int                            StartNode;
	int                            EndNode;
	TriMesh*                       themesh;
	map<int,int>                   Finish_id;
	vector<FuDijkstraNode>         FinishNode;
	priority_queue<FuDijkstraNode> WaitforSearch;
	void      RunDijkstra();
	void      output(vector<int>& os);
};


#endif