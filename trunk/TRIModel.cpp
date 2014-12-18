#include "TRIModel.h"

bool TRIModel::loadFromFile(const char* fileName){
	int aColorSet[] = {100, 100, 100, 100, 100, 100}; // front & back

	double max[3]={0.0, 0.0, 0.0};
	double min[3]={0.0, 0.0, 0.0};

	ifstream ifs(fileName);
	MyAssert(ifs.good());

	unsigned numF;
	ifs >> numF;

	triangleList.clear();
	for (unsigned f = 0; f < numF; f++) {
		Triangle tri(aColorSet);
		for (unsigned v = 0; v < 3; v++) {
			static DATA aData[6];
			ifs >> aData[0] >> aData[1] >> aData[2] 
				>> aData[3] >> aData[4] >> aData[5];
			
			for (unsigned i = 0; i < 3; i++) {
				if (aData[i] < min[i]) {
					min[i] = aData[i];
				} else {}
				if (aData[i] > max[i]) {
					max[i] = aData[i];
				} else {}
			}

			tri.loadVertex(v, aData);
		} // vertex
		triangleList.push_back(tri);
	} // face

	ifs.close();

	for(unsigned i = 0; i < 3; i++) {
		center[i] = (min[i] + max[i]) / 2;
	}
	return true;
}

/*
bool TRIModel::loadFromFile(const char* fileName){
	char tmp_string[100] = "";
	double max[3]={0.0, 0.0, 0.0};
	double min[3]={0.0, 0.0, 0.0};

	FILE* inFile = fopen(fileName, "r");
	if(!inFile){
		cout << "Can not open object File \"" << fileName << "\" !" << endl;
		return false;
	}

	cout <<"Loading \"" << fileName << "\" !" << endl;
	while(fscanf(inFile,"%s",tmp_string) != EOF){
		double tmp_double[6];
		int tmp_int[6];
		fscanf(inFile,"%d %d %d %d %d %d",&tmp_int[0],&tmp_int[1], &tmp_int[2], &tmp_int[3], &tmp_int[4], &tmp_int[5]);
		Triangle tmp_triangle(tmp_int);
		for(int i = 0; i < 3; i++){
			fscanf(inFile,"%lf %lf %lf %lf %lf %lf",&tmp_double[0],&tmp_double[1], &tmp_double[2], &tmp_double[3], &tmp_double[4], &tmp_double[5]);
			for(int j = 0; j < 3; j++){
				if(tmp_double[j] < min[j]){
					min[j] = tmp_double[j];
				}
				if(tmp_double[j] > max[j]){
					max[j] = tmp_double[j];
				}
			}
			tmp_triangle.loadVertex(i, tmp_double);
		}
		triangleList.push_back(tmp_triangle);
	}
	for(int i = 0; i < 3; i++){
		center[i] = (min[i] + max[i]) / 2;
	}
	return true;
}
*/
void TRIModel::copy(TRIModel * t){
	for(int i = 0; i < 3; i++){
		center[i] = t->center[i];
	}
	triangleList = t->triangleList;
}

TRIModel::TRIModel(){
}

TRIModel::~TRIModel(){
}