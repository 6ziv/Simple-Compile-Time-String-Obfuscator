# Simple Obfuscator

#### Introduction

a header-only compile-time AES-based string obfuscator.



I have been looking for a compile-time string obfuscator to use in my project. However, I cannot find one that still works with the latest compilers and meets all my needs. Thus I have decided to develop one on my own.

During development, I referred to SergeyBel's [AES implementation](https://github.com/SergeyBel/AES).

Random keys are used by default. 



You may use https://github.com/adamyaxley/Obfuscate as a simple xor obfuscator (which is far more efficient)



#### Usage

```C++
auto msg=Obfuscator::ProtectedArray<RANDOM_SEED,char,12>("Hello World.").decrypt();
std::cout<<msg.data()<<std::endl;
msg.clear();
```

or you may use it with other data types.

```
constexpr int raw_data[8]={0,1,2,3,4,5,6,7};
auto msg=Obfuscator::ProtectedArray<RANDOM_SEED,int,8>(raw_data).decrypt();
std::cout<<msg.data()<<std::endl;
msg.clear();
```

For convenience, some macros are provided. You can rewrite the examples as belows.

```
DECLARE_PROTECTED(char,12,encrypted_msg,"Hello World");
{
	std::cout<<encrypted_msg.decrypt().data()<<std::endl;
}
```

or

```
{
	std::cout<<PROTECTED(char,12,"Hello World")<<std::endl;
}
```

And, we have a convience macro defined for string literals.

```
std::cout<<PROTECTED_STRING("Hello World")<<std::endl;
```



Note that the decrypted resource will be cleared (filled with zero) before the temporary being destroyed, which means that the code section below

```
const char* msg=PROTECTED_STRING("Hello World");
std::cout<<msg<<std::endl;
```

will produce an empty output.



If you do not want these macros to ruin your code, you can define OBFUSCATOR_NO_MACROS to disable them, as well as the macro RANDOM_SEED. Then, you should manually specify a seed for each encrypted data.

Do not put two declarations in the same line. The macro RANDOM_SEED uses the `__LINE__` macro to generate randomness.



#### Requirements

C++20 is required to build this. To be accurate, I have been using 'consteval' here and there, and 'bit_cast' is also used. As for gcc, gcc 11 or higher is required.

When USE_BOOST is defined, boost will be required. The library will use boost's BOOST_FORCEINLINE on fetching and clearing decrypted data, as well as decrypting data, to make it harder to dump the decrypted data.



#### License

MIT