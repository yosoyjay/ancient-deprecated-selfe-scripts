function [db_conn] = connect_to_db(varargin)
% The function returns a handle to a database.  By default it connects
% to the CMOP database, but that can be changed via varargs.  The varargs
% are a skeleton of what may be required for other databases.
%
% Input: None if using CMOP database,   else
%   db_name   = varargin{1}
%   username = varargin{2}
%   password  = varargin{3}
%   db_url    = varargin{4}
%   db_driver = varargin{5}
%
%  lopezj - 7/2/2011
%

% Steps and settings to connect Matlab to the database
opt_args = size(varargin,2);
if opt_args == 5
	db_driver = varargin{5};
else
	db_driver = '/home/workspace/project/lopezj/bin/postgresql-9.0-801.jdbc4.jar';
end
javaaddpath(db_driver);

if opt_args == 4
	db_name = varargin{1};
	username = varargin{2};
	password = varargin{3};
	db_url = varargin{4};
else
	db_name   = 'cmop';
	username  = 'research';
	password  = 'rEsearch';
	db_driver = 'org.postgresql.Driver';
	db_url    = 'cdb02.stccmop.org:5432';
end

time_zone = '-08';
db_path  = sprintf('jdbc:postgresql://%s/%s', db_url, db_name);
db_conn  = database(db_name, username, password, db_driver, db_path);

