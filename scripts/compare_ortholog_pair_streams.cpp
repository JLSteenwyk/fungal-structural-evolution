// Exact merge of two sorted reciprocal-pair streams; see Python runner for provenance.
#include <array>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using Key = std::array<char,12>;
struct Stream {
    FILE* file;
    std::vector<char> buffer;
    size_t pos=0, size=0;
    Key previous{};
    bool have_previous=false;
    unsigned long long count=0;
    explicit Stream(const char* path):buffer(28*65536) {
        file=fopen(path,"rb"); if(!file) throw std::runtime_error("Cannot open stream");
    }
    ~Stream(){fclose(file);}
    bool next(Key& key) {
        if(pos==size){size=fread(buffer.data(),1,buffer.size(),file);pos=0;
            if(ferror(file))throw std::runtime_error("Read error");
            if(size%28)throw std::runtime_error("Truncated reciprocal pair");
            if(!size)return false;}
        const char* p=buffer.data()+pos;pos+=28;
        if(p[12]!='0'||p[13]!='\n'||p[26]!='1'||p[27]!='\n'||memcmp(p,p+14,12))
            throw std::runtime_error("Pair is not exactly once in each direction");
        for(int j=0;j<12;++j) if(!((p[j]>='0'&&p[j]<='9')||(p[j]>='a'&&p[j]<='f')))
            throw std::runtime_error("Invalid hexadecimal ID");
        if(memcmp(p,p+6,6)>=0)throw std::runtime_error("Noncanonical pair");
        memcpy(key.data(),p,12);
        if(have_previous && !(previous<key))throw std::runtime_error("Duplicate or unsorted pair");
        previous=key;have_previous=true;++count;return true;
    }
};
int main(int argc,char** argv){
    try {
        if(argc!=3)throw std::runtime_error("Usage: compare LEFT RIGHT");
        Stream a(argv[1]),b(argv[2]);Key x,y;bool ax=a.next(x),by=b.next(y);
        unsigned long long both=0,left=0,right=0;
        while(ax||by){
            if(ax&&by&&x==y){++both;ax=a.next(x);by=b.next(y);}
            else if(ax&&(!by||x<y)){++left;ax=a.next(x);}
            else {++right;by=b.next(y);}
        }
        std::cout<<"{\"left_pairs\":"<<a.count<<",\"right_pairs\":"<<b.count
                 <<",\"shared_pairs\":"<<both<<",\"left_only_pairs\":"<<left
                 <<",\"right_only_pairs\":"<<right<<",\"union_pairs\":"<<(both+left+right)
                 <<",\"symmetric_difference_pairs\":"<<(left+right)<<"}\n";
    } catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
