#include<bits/stdc++.h>
using namespace::std;

void parse(string s, string &a, string &b){
	a = "";
	b = "";
	bool toA = true;
	for(int i = 0; i < s.size(); i++){
		if(s[i] == '.'){
			toA = false;
			continue;
		}
		if(toA) a += s[i];
		else b += s[i];
	}
}

bool great(pair<string, string> a, pair<string, string> b){
	if(a.first.size() > b.first.size()) return true;
	if(a.first.size() < b.first.size()) return false;
	if(a.first > b.first) return true;
	if(a.first < b.first) return false;
	return a.second > b.second;
}

int main(){
	string re, im;
	pair<string, string> minre = make_pair("", ""), maxre = make_pair("", "");
	pair<string, string> minim = make_pair("", ""), maxim = make_pair("", "");
	while(cin >> re >> im){
		string reint, refloat;
		string imint, imfloat;
		parse(re, reint, refloat);
		parse(im, imint, imfloat);
		if(minre.first.empty() or great(minre, make_pair(reint, refloat))){
			minre = make_pair(reint, refloat);
		}
		if(maxre.first.empty() or great(make_pair(reint, refloat), maxre)){
			maxre = make_pair(reint, refloat);
		}
		if(minim.first.empty() or great(minim, make_pair(imint, imfloat))){
			minim = make_pair(imint, imfloat);
		}
		if(maxim.first.empty() or great(make_pair(imint, imfloat), maxim)){
			maxim = make_pair(imint, imfloat);
		}

	}
	cout << minre.first << "." << minre.second << ' ';
	cout << minim.first << "." << minim.second << endl;
	cout << maxre.first << "." << maxre.second << ' ';
	cout << maxim.first << "." << maxim.second << endl;
	return 0;
}
