////http://stackoverflow.com/questions/18939869/how-to-get-the-slope-of-a-linear-regression-line-using-c
////http://terpconnect.umd.edu/~toh/spectrum/CurveFitting.html
//http://www.yolinux.com/TUTORIALS/LinuxTutorialC++STL.html
//double slope(const std::vector<double>& x, const std::vector<double>& y) {
//    const double n    = x.size();
//    const double s_x  = std::accumulate(x.begin(), x.end(), 0.0);
//    const double s_y  = std::accumulate(y.begin(), y.end(), 0.0);
//    const double s_xx = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
//    const double s_xy = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
//    const double a    = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
//    return a;
//}


double slope(const ArrayType x, const ArrayType y) {

	const double n    = x.size();
	const double s_x  = std::accumulate(x.begin(), x.end(), 0.0);
	const double s_y  = std::accumulate(y.begin(), y.end(), 0.0);
	const double s_xx = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
	const double s_xy = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
	const double sl    = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);

	return sl;
}


double intercept(const ArrayType x, const ArrayType y) {

	const double n    = x.size();
	const double s_x  = std::accumulate(x.begin(), x.end(), 0.0);
	const double s_y  = std::accumulate(y.begin(), y.end(), 0.0);
	const double meanx =s_x/n;
	const double meany =s_y/n;
	const double s_xx = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
	const double s_xy = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
	const double sl    = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
	const double I = meany-(sl*meanx);
	return I;
}


double ssy(const ArrayType x, const ArrayType y) {

	const double n    = x.size();
	const double s_x  = std::accumulate(x.begin(), x.end(), 0.0);
	const double s_y  = std::accumulate(y.begin(), y.end(), 0.0);
	const double meanx =s_x/n;
	const double meany =s_y/n;
	const double s_xx = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
	const double s_xy = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
	const double sl    = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
	double yss;
	ArrayType::const_iterator s;
	for (s=y.begin();s!=y.end();s++)
		yss+=std::pow((*s)*meany,2);


	return yss;
}

double ssr(const ArrayType x, const ArrayType y) {

	const double n    = x.size();
	const double s_x  = std::accumulate(x.begin(), x.end(), 0.0);
	const double s_y  = std::accumulate(y.begin(), y.end(), 0.0);
	const double meanx =s_x/n;
	const double meany =s_y/n;
	const double s_xx = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
	const double s_xy = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
	const double sl    = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
	const double I = meany-(sl*meanx);
	ArrayType rss;
	//ArrayType::const_iterator s;
	for (int s=0;s<x.size();s++){
		rss.push_back(std::pow(y[s]-I-(x[s]*sl),2));
}
	double R=std::accumulate(rss.begin(),rss.end(),0.0);

	return R;
}


double r2(const ArrayType x, const ArrayType y) {

	const double n    = x.size();
	const double s_x  = std::accumulate(x.begin(), x.end(), 0.0);
	const double s_y  = std::accumulate(y.begin(), y.end(), 0.0);
	const double meanx =s_x/n;
	const double meany =s_y/n;
	const double s_xx = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
	const double s_xy = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
	const double sl    = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
	const double I = meany-(sl*meanx);
	ArrayType rss;
	//ArrayType::const_iterator s;
for (int s=0;s<x.size();s++){
		rss.push_back(std::pow(y[s]-I-(x[s]*sl),2));
}
	double R=std::accumulate(rss.begin(),rss.end(),0.0);
	double yss;
	ArrayType::const_iterator s;
	for (s=y.begin();s!=y.end();s++)
		yss+=std::pow((*s)*meany,2);


	return 1.0-(R/yss);
}



double sdslope(const ArrayType x, const ArrayType y) {

	const double n    = x.size();
	const double s_x  = std::accumulate(x.begin(), x.end(), 0.0);
	const double s_y  = std::accumulate(y.begin(), y.end(), 0.0);
	const double meanx =s_x/n;
	const double meany =s_y/n;
	const double s_xx = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
	const double s_xy = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
	const double sl    = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
	const double I = meany-(sl*meanx);
	ArrayType rss;
	//ArrayType::const_iterator s;
	for (int s=0;s<x.size();s++){
		rss.push_back(std::pow(y[s]-I-(x[s]*sl),2));
}
	double R=std::accumulate(rss.begin(),rss.end(),0.0);
double yss;
	ArrayType::const_iterator s;
	for (s=y.begin();s!=y.end();s++)
		yss+=std::pow((*s)*meany,2);


	return std::sqrt(R/(n-2))*std::sqrt(n/(n*s_xx-s_x*s_x));
}


double sdIntercept(const ArrayType x, const ArrayType y) {

	const double n    = x.size();
	const double s_x  = std::accumulate(x.begin(), x.end(), 0.0);
	const double s_y  = std::accumulate(y.begin(), y.end(), 0.0);
	const double meanx =s_x/n;
	const double meany =s_y/n;
	const double s_xx = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
	const double s_xy = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
	const double sl    = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
	const double I = meany-(sl*meanx);
	ArrayType rss;
	//ArrayType::const_iterator s;
	for (int s=0;s<x.size();s++){
		rss.push_back(std::pow(y[s]-I-(x[s]*sl),2));
}
	double R=std::accumulate(rss.begin(),rss.end(),0.0);
double yss;
	ArrayType::const_iterator s;
	for (s=y.begin();s!=y.end();s++)
		yss+=std::pow((*s)*meany,2);


	return std::sqrt(R/(n-2))*std::sqrt(s_xx/(n*s_xx-s_x*s_x));
}

