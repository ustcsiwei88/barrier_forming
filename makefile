SYS=LINUX
GUROBI_INCLUDE=/opt/gurobi910/linux64/include/
GUROBI_LIB=/opt/gurobi910/linux64/lib/

ifeq ($(SYS), WIN)
	CXX=cl
	CXXFLAGS=/permissive- /GS /GL /W3 /Gy /Zc:wchar_t /I"C:\gurobi900\win64\include" \
		/Zi /Gm- /O2 /Zc:inline /fp:precise /D "NDEBUG" \
		/D "_CONSOLE" /D "_MBCS" /errorReport:prompt /WX- /Zc:forScope /Gd /Oi /MT \
		/FC /EHsc /nologo  \
		/diagnostics:classic
	LINKFLAGS= /link "C:\gurobi900\win64\lib\gurobi_c++mt2017.lib" "C:\gurobi900\win64\lib\gurobi90.lib" \
                /LIBPATH:"C:\gurobi900\win64\lib"
else
	CXX=g++
	CXXFLAGS= -I$(GUROBI_INCLUDE) -Iinclude/ -O2 -m64 -std=c++17
	LINKFLAGS= -L$(GUROBI_LIB) -Wl,-rpath=$(GUROBI_LIB) -lgurobi_c++ -lgurobi91 -lm 
endif


all: min_seg min_seg_sample min_seg_multi_sets min_distance

min_seg: min_segment.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LINKFLAGS) 

min_seg_sample: min_segment_sample.cpp
	$(CXX) $(CXXFLAGS)  -o $@ $< $(LINKFLAGS) 

min_seg_multi_sets: min_segment_multi_sets.cpp
	$(CXX) $(CXXFLAGS)  -o $@ $< $(LINKFLAGS) 

min_distance: min_distance.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LINKFLAGS)



clean:
	rm -f min_distance min_seg_multi_sets min_seg_sample min_seg
