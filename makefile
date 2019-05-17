
SRCS := $(wildcard *.cpp)
BINS := $(SRCS:%.cpp=%)

.PHONY = all ${BINS}

test: ${BINS}
	@echo ${BINS}

test: ${BINS}
	@echo ${BINS}


%: %.o
	@echo "create"

%.o: %.cpp
	@echo "s"
#
# fd-plate: fd-plate.cpp
# 	clang++ -framework OpenAl fd-plate.cpp -o fd-plate.app
#
# hello:
# 	@echo "hello"
#
# clean:
# 	#rm .stdafx.c++.pch
# 	rm *.app
#
# file={1,2,3}
#
# list:
# 	for number in ${file} ; do \
#         echo $$number ; \
#   done


# RULE:	DEPENDENCY LINE
# [tab]ACTION LINE(S)

# DEPENDENCY LINE:					TARGET FILES:	SOURCE FILES
