extern void S_EMTP_WIN32();
extern void S_EMTP_linux();

#include <iostream>
using namespace std;
int main()
{
	//======linux Program======//
#ifdef linux
	S_EMTP_linux();
	cout<<"linux S_EMTP Program completed"<<endl;
#endif
	//======Win32  Program======//
#ifdef WIN32
	S_EMTP_WIN32();
	cout<<"WIN32 S_EMTP Program completed"<<endl;
#endif
	return 0;

}