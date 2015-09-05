#include "solid.h"
using namespace std;
class CFloyd
{
public:
	int INF;
	vector<vector<double> > floyd_Bigmap;
	vector<vector<double> > floyd_dist;
	vector<double> floyd_map;
	CFloyd(void);
	~CFloyd(void);
	void init(int n);
	void floyd(int n);
};