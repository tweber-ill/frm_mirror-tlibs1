find_path(Minuit2_INCLUDE_DIRS
	NAMES MnMigrad.h
	PATH_SUFFIXES root Minuit2 Minuit root/Minuit2 root/Minuit
	HINTS /usr/include/root/Minuit2 /usr/local/include/root/Minuit2 /opt/local/include/root/Minuit2
	DOC "Root/Minuit2 include directories"
)

# also include root base dir
list(APPEND Minuit2_INCLUDE_DIRS "${Minuit2_INCLUDE_DIRS}/..")


find_library(Minuit2_LIBRARIES
	NAMES Minuit2
	HINTS /usr/lib64/root /usr/lib/root /usr/lib32/root
	DOC "Minuit2 library"
)


message("Minuit include directories: ${Minuit2_INCLUDE_DIRS}")
message("Minuit library: ${Minuit2_LIBRARIES}")
