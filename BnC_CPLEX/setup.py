import os.path as opath
import platform
import os
import sys
from functools import reduce
from glob import glob

_platform = platform.system()

SYSTEM = 'x86-64_osx' if _platform == 'Darwin' else 'x86-64_linux'
DYL_EXT = 'dylib' if _platform == 'Darwin' else 'so'

CXX = 'clang++' if _platform == 'Darwin' else 'g++'
C = 'clang' if _platform == 'Darwin' else 'gcc'

CPLEX_LN_FLAGS_OSX = ['-lconcert', '-lilocplex', '-lcplex', 
						'-m64', '-lm', '-lpthread', 
						'-framework', 'CoreFoundation', 
						'-framework', 'IOKit',
						'-stdlib=libc++']
CPLEX_LN_FLAGS_LINUX = ['-lconcert', '-lilocplex', '-lcplex', 
						'-lm', '-lpthread', '-ldl']
CPLEX_LN_FLAGS = ' '.join(CPLEX_LN_FLAGS_OSX) \
					if _platform == 'Darwin' else \
				 ' '.join(CPLEX_LN_FLAGS_LINUX)

SRC_DIR = opath.join(os.getcwd(), 'src')
BUILD_DIR = opath.join(os.getcwd(), 'build')
LIB_DIR = opath.join(os.getcwd(), 'lib')

LIB_FORMAT = 'static_pic'
SRC_EXT, OBJ_EXT = 'cpp', 'o'

CXX_OPT = '-std=c++11 -DIL_STD'

UTIL_NAME = 'util'
###############################################################
CPLEX_HOME = None
CPLEX_DIR, CPLEX_INC_DIR, CPLEX_LIB_DIR = None, None, None
CONCERT_DIR, CONCERT_INC_DIR, CONCERT_LIB_DIR = None, None, None

UTIL_DIR = None
UTIL_INC_DIR, UTIL_LIB_DIR = None, None

###############################################################
MAIN_CPP = 'main.%s' % SRC_EXT
MAIN_OBJ = opath.join(BUILD_DIR, 'main.%s' % OBJ_EXT)

###############################################################
SOURCES = glob(SRC_DIR + '/**/*.%s' % SRC_EXT, recursive=True)
OBJECTS = [fpath.replace('src', 'build').replace(SRC_EXT, OBJ_EXT) for fpath in SOURCES]

###############################################################
SEQ_NAME='sequencer'
SEQ_DIR=opath.join(BUILD_DIR, 'GH')
SEQ_C=reduce(opath.join, [SRC_DIR, 'GH', '%s.c' % SEQ_NAME])
SEQ_OBJ=opath.join(SEQ_DIR, '%s.%s' % (SEQ_NAME, OBJ_EXT))
SEQ_STL=opath.join(LIB_DIR, 'lib%s.a' % SEQ_NAME)
SEQ_DYL=opath.join(LIB_DIR, 'lib%s.%s' % (SEQ_NAME, DYL_EXT))

###############################################################
CXX_INC_DIRS, CXX_LIB_DIRS = None, None
CXX_FLAGS, CXX_LN_FLAGS = None, None

###############################################################
PROGRAM = 'router'
PROGRAM_STL = opath.join(LIB_DIR, 'lib%s.a' % PROGRAM)

###############################################################

def hasFlag(flag):
	for i in range(len(sys.argv)):
		if sys.argv[i] == flag:
			return True
	else:
		return False

def valueOf(flag):
	for i in range(len(sys.argv)):
		if sys.argv[i] == flag:
			return sys.argv[i + 1]
	else:
		return None

def set_paths():
	global CPLEX_HOME
	global CPLEX_DIR, CPLEX_INC_DIR, CPLEX_LIB_DIR
	global CONCERT_DIR, CONCERT_INC_DIR, CONCERT_LIB_DIR
	if hasFlag('-cp'):
		CPLEX_HOME = valueOf('-cp')
	else:
		if platform.system() == 'Darwin':
			CPLEX_HOME = '/Applications/CPLEX_Studio128'
		else:
			print('Please a path directing CPLEX_HOME')
			assert False
	CPLEX_DIR = opath.join(CPLEX_HOME, 'cplex')
	CPLEX_INC_DIR = opath.join(CPLEX_DIR, 'include')
	CPLEX_LIB_DIR = reduce(opath.join, 
							[CPLEX_DIR, 'lib', SYSTEM, LIB_FORMAT])
	CONCERT_DIR = opath.join(CPLEX_HOME, 'concert')
	CONCERT_INC_DIR = opath.join(CONCERT_DIR, 'include')
	CONCERT_LIB_DIR = reduce(opath.join, 
							[CONCERT_DIR, 'lib', SYSTEM, LIB_FORMAT])
	#
	global UTIL_DIR
	global UTIL_INC_DIR, UTIL_LIB_DIR
	if hasFlag('-ut'):
		UTIL_DIR = valueOf('-ut')
	else:
		if platform.system() == 'Darwin':
			UTIL_DIR = '/Users/ckhan/workspace/util'
		else:
			print('Please a path directing UTIL_DIR')
			assert False
	UTIL_INC_DIR = opath.join(UTIL_DIR, 'include')
	UTIL_LIB_DIR = opath.join(UTIL_DIR, 'lib')
	#

def set_flags():
	if CPLEX_HOME is None:
		set_paths()
	#
	global CXX_INC_DIRS, CXX_LIB_DIRS
	global CXX_FLAGS, CXX_LN_FLAGS
	CXX_INC_DIRS = ' '.join(['-I%s' % dpath for dpath in 
				[CPLEX_INC_DIR, CONCERT_INC_DIR, UTIL_INC_DIR]])
	CXX_LIB_DIRS = ' '.join(['-L%s' % dpath for dpath in 
				[CPLEX_LIB_DIR, CONCERT_LIB_DIR, UTIL_LIB_DIR, LIB_DIR]])
	CXX_FLAGS = '%s %s' % (CXX_OPT, CXX_INC_DIRS)
	CXX_LN_FLAGS = '%s -l%s -l%s' % (CPLEX_LN_FLAGS, SEQ_NAME, UTIL_NAME)

def comCPP(src_fpath, obj_fpath):
	if not opath.exists(obj_fpath) or \
		opath.getctime(obj_fpath) < opath.getmtime(src_fpath):
		os.makedirs(opath.dirname(obj_fpath), exist_ok=True)
		print("Start compile") 
		print("Target file:", obj_fpath, "; Pre-req:", src_fpath) 
		res = os.system('%s -c %s %s -o %s' % (CXX, CXX_FLAGS, src_fpath, obj_fpath))
		assert res == 0
		print("Build Success!!!\n")

def comBASE():
	comCPP(MAIN_CPP, MAIN_OBJ)
	for i in range(len(SOURCES)):
		src_fpath, obj_fpath = SOURCES[i], OBJECTS[i]
		comCPP(src_fpath, obj_fpath)

def genSTL_ROUTER():
	comBASE()
	#
	if not opath.exists(PROGRAM_STL):
		os.makedirs(opath.dirname(PROGRAM_STL), exist_ok=True)
		print("Build a router static library")
		res = os.system("ar rvs %s %s" \
						% (PROGRAM_STL, ' '.join(OBJECTS)))
		assert res == 0
		print("Successfully build the library\n") 
	else:
		print("The library exists already\n") 

def genSTL_SEQ():
	if not opath.exists(SEQ_STL):
		os.makedirs(opath.dirname(SEQ_STL), exist_ok=True)
		print("Build a sequencer static library") 
		res = os.system("%s -c %s -o %s" \
						% (C, SEQ_C, SEQ_OBJ))
		assert res == 0
		res = os.system("ar rvs %s %s" \
						% (SEQ_STL, SEQ_OBJ))
		assert res == 0
		print("Successfully build the library\n") 
	else:
		print("The library exists already\n") 

def genBIN():
	comBASE()
	genSTL_SEQ()
	#
	cmd = "%s %s %s %s %s %s -o %s" \
					% (CXX, CXX_FLAGS, CXX_LIB_DIRS, CXX_LN_FLAGS,
						MAIN_OBJ, ' '.join(OBJECTS), PROGRAM)
	res = os.system(cmd)
	assert res == 0
	

if __name__ == '__main__':
	set_paths()
	set_flags()
	if len(sys.argv) == 1:
		genBIN()
		genSTL_ROUTER()
	else:
		if hasFlag('clean'):
			os.system('rm -rf %s' % BUILD_DIR)
			os.system('rm -rf %s' % LIB_DIR)
			os.system('rm %s' % PROGRAM)
		if hasFlag('stl'):
			genSTL_ROUTER()
			genSTL_SEQ()



	

