#include <stdio.h>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <fstream>
#include <typeinfo>

using namespace std;

//class Convert{
//public:
//	// Convert T, which should be a primitive, to a string.
//	template <typename T>
//	static string T_to_string(T const &val) {
//		ostringstream ostr;
//		ostr << val;
//
//		return ostr.str();
//	}
	//
//	// Convert a string to T.	
//	template <typename T>
//	static T string_to_T(string const &val) {
//		istringstream istr(val);
//		T returnVal;
//		if (!(istr >> returnVal))
//			//exitWithError("CFG: Not a valid " + (string)typeid(T).name() + " received!\n");
//
//		return returnVal;
//	}
//
//	template <typename T>
//	static string string_to_T(string const &val){
//		return val;
//	}
//};

class ConfigFile{
private:
	map<string, string> contents;
	string fName;
	void removeComment(string &line) const{
    	if (line.find(';') != line.npos)
		line.erase(line.find(';'));
	}
	bool onlyWhitespace(const string &line) const{
    	return (line.find_first_not_of(' ') == line.npos);
	}
	bool validLine(const string &line) const{
		string temp = line;
		temp.erase(0, temp.find_first_not_of("\t "));
		if (temp[0] == '='){
			return false;
		}
		for (size_t i = temp.find('=') + 1; i < temp.length(); i++){	
			if (temp[i] != ' '){
				return true;
			}
		}
		return false;
	}
	void extractKey(string &key, size_t const &sepPos, const string &line) const{
    	key = line.substr(0, sepPos);
    	if (key.find('\t') != line.npos || key.find(' ') != line.npos){
			key.erase(key.find_first_of("\t "));
		}
	}
	void extractValue(string &value, size_t const &sepPos, const string &line) const{
		value = line.substr(sepPos + 1);
		value.erase(0, value.find_first_not_of("\t "));
		value.erase(value.find_last_not_of("\t ") + 1);
	}
	void extractContents(const string &line) {
		string temp = line;
	        // Erase leading whitespace from the line.
		temp.erase(0, temp.find_first_not_of("\t "));
		size_t sepPos = temp.find('=');

		string key, value;
		extractKey(key, sepPos, temp);
		extractValue(value, sepPos, temp);

		if (!keyExists(key))
			contents.insert(pair<string, string>(key, value));
		//else
			//exitWithError("CFG: Can only have unique key names!\n");
	}

	// lineNo = the current line number in the file.
	// line = the current line, with comments removed.
	void parseLine(const string &line, size_t const lineNo){
		if (line.find('=') == line.npos){
			//exitWithError("CFG: Couldn't find separator on line: " + Convert::T_to_string(lineNo) + "\n");
		}
		if (!validLine(line)){
			//exitWithError("CFG: Bad format for line: " + Convert::T_to_string(lineNo) + "\n");
		}
			extractContents(line);
	}
	void ExtractKeys(){
		ifstream file;
		file.open(fName.c_str());
		if (!file){
			//exitWithError("CFG: File " + fName + " couldn't be found!\n");
		}
		string line;
		size_t lineNo = 0;
		while (getline(file, line)){
			lineNo++;
			string temp = line;
			if (temp.empty()){
				continue;
			}
			removeComment(temp);
			if (onlyWhitespace(temp)){
				continue;
			}
			parseLine(temp, lineNo);
		}
		file.close();
	}

public:
	ConfigFile(const string &fName){
		this->fName = fName;
		ExtractKeys();
	}
	bool keyExists(const string &key) const{
		return contents.find(key) != contents.end();
	}

	template <typename ValueType>
	ValueType getValueOfKey(const string &key, ValueType const &defaultValue = ValueType()) const{
		if (!keyExists(key)){
			return defaultValue;
		}
		//return Convert::string_to_T<ValueType>(contents.find(key)->second);
		return contents.find(key)->second;
	}
};

bool bool_convert(string option){
	if(option == "true"){
		return true;
	} else if(option == "false"){
		return false;
	} else{
		return 0;
	}
}

/* This function can be used to read a data file and plot the data contained. */

int main(){
	ConfigFile cfg("graph_config.cfg");
	
	string tmp;

	string input_filename_default = cfg.getValueOfKey<string>("input_filename_default");
	
	tmp = cfg.getValueOfKey<string>("output_to_latex");
	bool output_to_latex = bool_convert(tmp);
	
	tmp = cfg.getValueOfKey<string>("output_to_pdf");
	bool output_to_pdf = bool_convert(tmp);
		tmp = cfg.getValueOfKey<string>("output_filename_default");
		const char * output_filename_default = tmp.c_str();
	
	tmp = cfg.getValueOfKey<string>("use_auto_ranges");
	bool use_auto_ranges = bool_convert(tmp);
		tmp = cfg.getValueOfKey<string>("xrange");
		const char * xrange = tmp.c_str();
		tmp = cfg.getValueOfKey<string>("yrange");
		const char * yrange = tmp.c_str();

	tmp = cfg.getValueOfKey<string>("join_points_with_lines");
	bool join_points_with_lines = bool_convert(tmp);
	
	tmp = cfg.getValueOfKey<string>("draw_lines_and_points");
	bool draw_lines_and_points = bool_convert(tmp);
	
	tmp = cfg.getValueOfKey<string>("draw_legend");
	bool draw_legend = bool_convert(tmp);
	
	tmp = cfg.getValueOfKey<string>("title");
	const char * title = tmp.c_str();
		tmp = cfg.getValueOfKey<string>("xlabel");
		const char * xlabel = tmp.c_str();
		tmp = cfg.getValueOfKey<string>("ylabel");
		const char * ylabel = tmp.c_str();

	tmp = cfg.getValueOfKey<string>("apply_fit");
	bool apply_fit = bool_convert(tmp);
		tmp = cfg.getValueOfKey<string>("fit_equation");
		const char * fit_equation = tmp.c_str();
		tmp = cfg.getValueOfKey<string>("provide_guesses");
		bool provide_guesses = bool_convert(tmp);
			tmp = cfg.getValueOfKey<string>("variables");
			const char * variables = tmp.c_str();
			tmp = cfg.getValueOfKey<string>("initial_m");
			const char * initial_m = tmp.c_str();
			tmp = cfg.getValueOfKey<string>("initial_c");
			const char * initial_c = tmp.c_str();

	string input_filename;
	string output_filename;

	cout << "Enter Input filename: ";
	cin >> input_filename;
	if(input_filename=="1"){
		input_filename = input_filename_default;
		cout << input_filename << endl;
	}
	cout << "Enter Output filename: ";
	cin >> output_filename;
	if(output_filename=="1"){
		output_filename = output_filename_default;
		cout << output_filename << endl;
	}

	const char * with_lines = "";
	if(join_points_with_lines){
		with_lines = "w l ";
	}
	string lines_and_points = "";
	if(draw_lines_and_points){
		lines_and_points = ", \"" + input_filename + "\"";
	}

	const char *lines_and_points_char = lines_and_points.c_str();
	const char *input_filename_char = input_filename.c_str();
	const char *output_filename_char = output_filename.c_str();

	FILE *gnuplotPipe = popen("gnuplot -persist","w");

	// If gnuplot is found
	if (gnuplotPipe) {

		fprintf(gnuplotPipe, "%s%s%s","set title \"",title,"\"\n");
		fprintf(gnuplotPipe, "%s%s%s","set xlabel \"",xlabel,"\"\n");
		fprintf(gnuplotPipe, "%s%s%s","set ylabel \"",ylabel,"\"\n");
		
		if(output_to_latex){
			fprintf(gnuplotPipe, "%s","set terminal latex\n");
			fprintf(gnuplotPipe, "%s%s%s","set output \"", output_filename_char,".tex\"\n");
		}
		
		if(draw_legend){
			fprintf(gnuplotPipe,"set key default \n");
		}else{
			fprintf(gnuplotPipe,"unset key  \n");
		}

 		if(!use_auto_ranges){
 			fprintf(gnuplotPipe, "%s%s%s%s%s", "set xrange[",xrange,"] \nset yrange[",yrange,"] \n");
 		}

 		if(apply_fit){
			fprintf(gnuplotPipe, "%s%s", fit_equation,"\n");
			if(provide_guesses){
				fprintf(gnuplotPipe, "%s%s%s%s", initial_m,"\n", initial_c,"\n");
			}
			fprintf(gnuplotPipe, "%s%s%s%s%s", "fit f(x) \"",input_filename_char,"\" via ",variables,"\n");
 		}

		fprintf(gnuplotPipe, "%s%s%s%s%s","plot \"",input_filename_char,"\"", with_lines, lines_and_points_char);

		// fprintf(gnuplotPipe, "%s%d%s\n","plot \"data",x,".txt\"");

		// if(join_points_with_lines){
		// 	if(draw_lines_and_points){
		// 		fprintf(gnuplotPipe, "%s%s%s",", \"",filename,"\" w l");
		// 	}
		// }else{
		// 	fprintf(gnuplotPipe, "%s\n"," ");
		// }

		if(apply_fit){
			fprintf(gnuplotPipe, "%s",", f(x)");
		}
		fprintf(gnuplotPipe, "%s","\n");

		if(output_to_pdf){
			fprintf(gnuplotPipe, "%s","set term post enhanced color solid \"Helvetica\" 16\n");
			fprintf(gnuplotPipe, "%s%s%s","set output '| epstopdf --filter --outfile=",output_filename_char,".pdf'\n");
			fprintf(gnuplotPipe, "%s","replot\n");
		}

		fflush(stderr);

		// exit gnuplot
		fprintf(gnuplotPipe,"exit \n");
		pclose(gnuplotPipe);

	//If gnuplot is not installed
	}else{
		cout << endl << "ERROR: Could not open GnuPlot" << endl;
	}
}
