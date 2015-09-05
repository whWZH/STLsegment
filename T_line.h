class T_line  
{
public:
	long p0;
	long p1;
	int useCount;
	bool	friend operator==(T_line line0,T_line line1);

	T_line();
	virtual ~T_line();

};