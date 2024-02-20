#include <algorithm> // for copy
#include <iterator> // for ostream_iterator
#include <vector>
#include<fstream>
#include<sstream>
#include<iostream>
#include<ostream>
#include <string>


char no='N';

void bar(){
	std::cout<<"\n------------------------------------------"<<std::endl;
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
	}};


void firma(){
	std::cout << "\n\n\neros.montin@gmail.com" << std::endl;
	bar();

};

void checkIn(std::vector<std::string> & P,char *argv[],int argc){
	//std::vector<std::string>p;
	bar();
	std::cout<<"Parameters input"<<std::endl;
	bar();
	for (int o=1; o<argc;o++)
		std::cout<< "Param "<<o<<":\t"<< P[o]<<":\t"<<argv[o]<<std::endl;
	bar();
};


void readParametersFromStdin(int * P, int N, char *argv[],int argc, char no ){
	if (argc>N){
		if( *argv[N]!=no){*P=atoi(argv[N]);}}
}


void readParametersFromStdin(double * P, int N, char *argv[],int argc, char no ){
	if (argc>N){
		if( *argv[N]!=no){*P=atof(argv[N]);}}
};

void readParametersFromStdin(float * P, int N, char *argv[],int argc, char no ){
	if (argc>N){
		if( *argv[N]!=no){*P=atof(argv[N]);}}
};


void readParametersFromStdin(bool * P, int N, char *argv[],int argc, char no ){
	if (argc>N){
		if( *argv[N]!=no){
			if (atoi(argv[N])<1){
				*P=false;
				}else{
					*P=true;
					}
			};
		};
};


void readParametersFromStdin(std::string * P, int N, char *argv[],int argc, char no ){
	if (argc>N){
		if( *argv[N]!=no){*P=argv[N];}}
}

void pMetric(std::string VV,std::vector<std::string> & P,std::vector<std::string> & V){
	bar();
	std::cout<<VV<<" parameters "<<std::endl;
	bar();

	for (long unsigned int o=0; o<P.size();o++)
		std::cout<< P[o]<<":\t"<< V[o]<<std::endl;
	bar();
};

template <typename T> 
double toRad(const T& t){
double r=vnl_math::pi/180;
return (double)t*r;
};

template<typename T>
void printVector(std::vector<T> & v){
	std::copy(v.begin(),
              v.end(),
              std::ostream_iterator<T>(std::cout,",")
            );
};


template<typename T>
void printVectorLine(std::vector<T> & vect){
	typename std::vector<T>::iterator it;

	for (it= vect.begin(); it != vect.end(); it++){
		std::cout<<*it<<",";
	};
	std::cout<<"\n";
}

template<typename T>
void printVectorLine(std::string W, std::vector<T> & vect,std::string F ,double O){
	typename std::vector<T>::iterator it;
	std::cout<< W <<": ";
	for (it= vect.begin(); it != vect.end(); it++){
		std::cout<<*it<<", ";
	};

	std::cout<<F<<": "<<O<<std::endl;
}


template<typename T>
void printVectorToFile(std::vector<T> v, std::string txt){
std::ofstream output_file(txt.c_str());
 std::ostream_iterator<T> output_iterator(output_file, "\n");
 std::copy(v.begin(), v.end(), output_iterator);
}
