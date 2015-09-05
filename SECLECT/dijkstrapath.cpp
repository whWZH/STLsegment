#include "dijkstrapath.h"

typedef map<int,int>::iterator miter;

void DijkstraPath::RunDijkstra()
{
	//if uninitialized 
	if (!themesh || StartNode == -1 || EndNode == -1)
		return;

	FuDijkstraNode NodeS;
	NodeS.id = StartNode; NodeS.dist = 0;
	WaitforSearch.push(NodeS);

	while (!WaitforSearch.empty())
	{
		FuDijkstraNode tmpNode = WaitforSearch.top();
		WaitforSearch.pop();
	    
		pair<miter,bool> notInFinish;
		int fnode = FinishNode.size();      //index of the node 
		notInFinish = Finish_id.insert(make_pair(tmpNode.id,fnode));
		if (notInFinish.second)
		{
			FinishNode.push_back(tmpNode);
		}
		else
		{
			continue;
		}

		if (tmpNode.id == EndNode)
			break;

		vector<int>& neighbor = themesh->neighbors[tmpNode.id];
		for (unsigned i = 0; i < neighbor.size(); i++)
		{
			if (Finish_id.count(neighbor[i]))
				continue;
			FuDijkstraNode newNode;  
			newNode.id  = neighbor[i];
			newNode.pre = fnode;
            vec length  = themesh->vertices[newNode.id] - themesh->vertices[tmpNode.id];
			newNode.dist  = tmpNode.dist + len2(length);
			WaitforSearch.push(newNode);
		}
	}

}

void DijkstraPath::output(vector<int>& os)
{
	int end_id = -1;
	map<int,int>::iterator it = Finish_id.find(EndNode);
	if (it != Finish_id.end())
	{
		end_id = it->second;
	}
	else
		return;

	int pre_id = FinishNode[end_id].pre;
	if(pre_id==-1)return;
	while (FinishNode[pre_id].id != StartNode)
	{
		os.push_back(FinishNode[pre_id].id);
		pre_id = FinishNode[pre_id].pre;
	} 

}