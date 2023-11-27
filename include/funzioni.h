
#include <algorithm> // for copy
#include <iterator> // for ostream_iterator
#include <vector>
#include<fstream>
#include<sstream>
#include<iostream>



char no='N';

void bar(){
	std::cout<<"------------------------------------------"<<std::endl;
}
template <typename T> std::string tostr(const T& t) {
	std::ostringstream os;
	os<<t;
	return os.str();
}

void printVersion(){
	std::cerr<<"Verion : " << VERSION_MAJOR <<"." << VERSION_MINOR<< "\n"<<std::endl;
}



void printParam(std::vector<std::string> & A){
	bar();
	//	std::copy(command.begin(), command.end(), std::ostream_iterator<std::string>(std::cout, "\n "));
	std::vector<std::string>::iterator it;
	int a=-1;
	for (it=A.begin();it!=A.end();it++){
		a=a+1;
		std::cout<< "Param "<< a <<":\t"<< *it<<"\n";
	}}


void firma(){
	std::cout << "\n\n\neros.montin@polimi.it" << std::endl;
	bar();

}

void checkIn(std::vector<std::string> & P,char *argv[],int argc){
	//std::vector<std::string>p;
	bar();
	std::cout<<"Parameters input"<<std::cout;
	bar();
	for (int o=1; o<argc;o++)
		std::cout<< "Param "<<o<<":\t"<< P[o]<<":\t"<<argv[o]<<std::endl;
	bar();
}


void readParametersFromStdin(int * P, int N, char *argv[],int argc, char no ){
	if (argc>N){
		if( *argv[N]!=no){*P=atoi(argv[N]);}}
}


void readParametersFromStdin(double * P, int N, char *argv[],int argc, char no ){
	if (argc>N){
		if( *argv[N]!=no){*P=atof(argv[N]);}}
}

void readParametersFromStdin(float * P, int N, char *argv[],int argc, char no ){
	if (argc>N){
		if( *argv[N]!=no){*P=atof(argv[N]);}}
}





void pMetric(std::string VV,std::vector<std::string> & P,std::vector<std::string> & V){
	bar();
	std::cout<<VV<<" parameters "<<std::endl;
	bar();

	for (int o=0; o<P.size();o++)
		std::cout<< P[o]<<":\t"<< V[o]<<std::endl;
	bar();
}

template <typename T> 
double toRad(const T& t){
double r=vnl_math::pi/180;
return (double)t*r;
};



double signal2T1(double si, double so, double tr, double t10, double fa ){
double D,E10,A,B;
E10=vcl_exp(-tr/t10);
B=(1-E10)/(1-vcl_cos(toRad(fa))*E10);
A=B*si/so;
return 1/((-1/tr)* vcl_log((1-A)/(1-(vcl_cos(toRad(fa))*A))));
};

double T1t2GDc(double t1t, double t10, double re){
return ((1/t1t)-(1/t10))/re;
};




// reading a text file
//#include <iostream>
//#include <fstream>
//#include <string>
//using namespace std;



/*
std::vector<std::vector<double> > numbers;

std::string temp;

while (std::getline(infile, temp)) {
    std::istringstream buffer(temp);
    std::vector<double> line((std::istream_iterator<double>(buffer)),
                             std::istream_iterator<double>());

    numbers.push_back(line);
}

*/
