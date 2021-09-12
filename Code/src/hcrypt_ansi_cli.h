#ifndef __ANSI_CLI__
#define __ANSI_CLI__

#include <vector>

using namespace std;

class t_string;
class t_config;

/*
Portable command line interface. Defines universal usage of command line interface options.
*/
class t_ansi_cli
{
public:
	// Option prefix allows for customizing the prefix that specifies the options. For instance it can be "--" or "=" or "-".
	t_ansi_cli(int argc, char* argv[], const char* _option_prefix);
	~t_ansi_cli();

	char* get_value_by_option(const char* option, bool& success);
	bool is_flag_set(char* flag);
	vector<char*>* get_input_list();

	char* option_prefix;

private:
	char* exe_cmd;

	bool is_option(char* arg);

	// Options and corresponding values. Everything in the command line 
	vector<char*>* options;
	vector<char*>* values;	

	// Input: These are the values that are input to the program that do not have any options in front of them.
	// For example file names.
	vector<char*>* inputs;
};

bool get_cli_value(t_config* config, t_ansi_cli* cli, const char* config_id, const char* cli_id, char* val_buff, bool config_mandatory);


#endif // __ANSI_CLI__