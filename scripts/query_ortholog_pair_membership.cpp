// Complete reciprocal-stream scan for exact membership of sorted query pairs.
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
        if(argc!=4)throw std::runtime_error("Usage: membership STREAM QUERIES OUTPUT");
        Stream stream(argv[1]);
        FILE* queries=fopen(argv[2],"rb"); FILE* output=fopen(argv[3],"wb");
        if(!queries||!output)throw std::runtime_error("Cannot open query/output");
        Key target{}, previous{}; bool have_previous=false; unsigned long long requested=0,found=0;
        auto query_next=[&](){
            char line[13]; size_t n=fread(line,1,13,queries);
            if(ferror(queries))throw std::runtime_error("Query read error");
            if(!n)return false;
            if(n!=13||line[12]!='\n')throw std::runtime_error("Malformed query");
            for(int j=0;j<12;++j)if(!((line[j]>='0'&&line[j]<='9')||(line[j]>='a'&&line[j]<='f')))throw std::runtime_error("Invalid query ID");
            if(memcmp(line,line+6,6)>=0)throw std::runtime_error("Noncanonical query");
            memcpy(target.data(),line,12);
            if(have_previous&&!(previous<target))throw std::runtime_error("Queries not unique and sorted");
            previous=target;have_previous=true;++requested;return true;
        };
        auto write=[&](bool present){
            if(fwrite(target.data(),1,12,output)!=12||fprintf(output,"\t%d\n",int(present))<0)throw std::runtime_error("Write error");
            found+=present;
        };
        bool pending=query_next();Key current;
        while(stream.next(current)){
            while(pending&&target<current){write(false);pending=query_next();}
            if(pending&&target==current){write(true);pending=query_next();}
        }
        while(pending){write(false);pending=query_next();}
        if(fclose(output)!=0)throw std::runtime_error("Close error");fclose(queries);
        std::cout<<"{\"stream_pairs\":"<<stream.count<<",\"queries\":"<<requested<<",\"present\":"<<found<<",\"absent\":"<<requested-found<<"}\n";
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
