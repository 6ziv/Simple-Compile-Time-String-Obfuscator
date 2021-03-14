#include "SimpleObfuscator.hpp"
#include <iostream>
int main(){
	std::cout<<PROTECTED(char,12,"Hello World")<<std::endl;
	std::cout<<PROTECTED_STRING("Hello World2")<<std::endl;
	DECLARE_PROTECTED(char,100,str,"This is a header-only compile-time AES obfuscator. In order to use this library, c++20 is required.")
	auto decrypted = str.decrypt();
	std::cout<<decrypted.data()<<std::endl;
	decrypted.clear();
	
	auto msg2=SimpleObfuscator::ProtectedArray<RANDOM_SEED,char,109>("Though it will work without external requirements, it is recommended that users use this library with boost.").decrypt();
	std::cout<<msg2.data()<<std::endl;
	return 0;
}
