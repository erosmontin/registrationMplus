
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


std::string availableMetrics(){
std::string S=	"metric mode (0 MATTES, 1 NGF, 2 MANGF+, 3 MANG*, 4 MANGFMSE, MSE)";
	return S;
}

