#pragma once
#include <array>
#include <bit>
#include <cstdint>
#include <cstring>
#include <random>
#ifdef USE_BOOST
#include <boost/config/detail/suffix.hpp>
#define OBFUSCATOR_INLINE BOOST_FORCEINLINE
#else
#define OBFUSCATOR_INLINE inline
#endif
namespace SimpleObfuscator {
	namespace Impl {
		constexpr uint8_t SBox[16][16] = {
			0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5,
			0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,
			0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0,
			0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,
			0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc,
			0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
			0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a,
			0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,
			0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0,
			0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,
			0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b,
			0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
			0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85,
			0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,
			0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5,
			0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,
			0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17,
			0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
			0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88,
			0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,
			0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c,
			0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,
			0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9,
			0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
			0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6,
			0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,
			0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e,
			0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,
			0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94,
			0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
			0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68,
			0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16
		};

		constexpr uint8_t invSBox[16][16] = {
			0x52, 0x09, 0x6a, 0xd5, 0x30, 0x36, 0xa5, 0x38,
			0xbf, 0x40, 0xa3, 0x9e, 0x81, 0xf3, 0xd7, 0xfb,
			0x7c, 0xe3, 0x39, 0x82, 0x9b, 0x2f, 0xff, 0x87,
			0x34, 0x8e, 0x43, 0x44, 0xc4, 0xde, 0xe9, 0xcb,
			0x54, 0x7b, 0x94, 0x32, 0xa6, 0xc2, 0x23, 0x3d,
			0xee, 0x4c, 0x95, 0x0b, 0x42, 0xfa, 0xc3, 0x4e,
			0x08, 0x2e, 0xa1, 0x66, 0x28, 0xd9, 0x24, 0xb2,
			0x76, 0x5b, 0xa2, 0x49, 0x6d, 0x8b, 0xd1, 0x25,
			0x72, 0xf8, 0xf6, 0x64, 0x86, 0x68, 0x98, 0x16,
			0xd4, 0xa4, 0x5c, 0xcc, 0x5d, 0x65, 0xb6, 0x92,
			0x6c, 0x70, 0x48, 0x50, 0xfd, 0xed, 0xb9, 0xda,
			0x5e, 0x15, 0x46, 0x57, 0xa7, 0x8d, 0x9d, 0x84,
			0x90, 0xd8, 0xab, 0x00, 0x8c, 0xbc, 0xd3, 0x0a,
			0xf7, 0xe4, 0x58, 0x05, 0xb8, 0xb3, 0x45, 0x06,
			0xd0, 0x2c, 0x1e, 0x8f, 0xca, 0x3f, 0x0f, 0x02,
			0xc1, 0xaf, 0xbd, 0x03, 0x01, 0x13, 0x8a, 0x6b,
			0x3a, 0x91, 0x11, 0x41, 0x4f, 0x67, 0xdc, 0xea,
			0x97, 0xf2, 0xcf, 0xce, 0xf0, 0xb4, 0xe6, 0x73,
			0x96, 0xac, 0x74, 0x22, 0xe7, 0xad, 0x35, 0x85,
			0xe2, 0xf9, 0x37, 0xe8, 0x1c, 0x75, 0xdf, 0x6e,
			0x47, 0xf1, 0x1a, 0x71, 0x1d, 0x29, 0xc5, 0x89,
			0x6f, 0xb7, 0x62, 0x0e, 0xaa, 0x18, 0xbe, 0x1b,
			0xfc, 0x56, 0x3e, 0x4b, 0xc6, 0xd2, 0x79, 0x20,
			0x9a, 0xdb, 0xc0, 0xfe, 0x78, 0xcd, 0x5a, 0xf4,
			0x1f, 0xdd, 0xa8, 0x33, 0x88, 0x07, 0xc7, 0x31,
			0xb1, 0x12, 0x10, 0x59, 0x27, 0x80, 0xec, 0x5f,
			0x60, 0x51, 0x7f, 0xa9, 0x19, 0xb5, 0x4a, 0x0d,
			0x2d, 0xe5, 0x7a, 0x9f, 0x93, 0xc9, 0x9c, 0xef,
			0xa0, 0xe0, 0x3b, 0x4d, 0xae, 0x2a, 0xf5, 0xb0,
			0xc8, 0xeb, 0xbb, 0x3c, 0x83, 0x53, 0x99, 0x61,
			0x17, 0x2b, 0x04, 0x7e, 0xba, 0x77, 0xd6, 0x26,
			0xe1, 0x69, 0x14, 0x63, 0x55, 0x21, 0x0c, 0x7d
		};
		enum { BlockSize = 16 };
		typedef std::array<uint8_t, BlockSize> Block;
		template <size_t N>
		inline consteval Block padding(std::array<uint8_t, N> x) {
			Block result{};
			std::size_t index = 0;
			for (auto& item : x) {
				result[index] = std::move(item);
				++index;
			}
			while (index < BlockSize) {
				result[index] = 0;
				++index;
			}
			return result;
		}
		inline constexpr Block xorBlock(Block x1, Block x2)
		{
			Block result{};
			for (std::size_t i = 0; i < BlockSize; i++)
			{
				result[i] = static_cast<uint8_t>(x1[i] ^ x2[i]);
			}
			return result;
		}
		inline void xorBlockInplace(uint8_t* x1, Block x2)
		{
			for (std::size_t i = 0; i < BlockSize; i++)
				x1[i] ^= x2[i];
		}
		template <size_t keyLength>
		class parameters {
		public:
			enum { Nk = 0 };
			enum { Nr = 0 };
		};
		template<>
		class parameters<128> {
		public:
			enum { Nk = 4 };
			enum { Nr = 10 };
		};
		template<>
		class parameters<192> {
		public:
			enum { Nk = 6 };
			enum { Nr = 12 };
		};
		template<>
		class parameters<256> {
		public:
			enum { Nk = 8 };
			enum { Nr = 14 };
		};
		using StatePart = std::array<uint8_t, 4>;
		using State = std::array<StatePart, 4>;
		using RoundKeyPart = Block;

		template <size_t keyLength>
		using RoundKey = std::array<RoundKeyPart, parameters<keyLength>::Nr + 1>;

		inline consteval State AddRoundKey(State state, RoundKeyPart key)
		{
			State result{};
			for (size_t i = 0; i < 4; i++) {
				for (size_t j = 0; j < 4; j++)
				{
					result[i][j] = static_cast<uint8_t>(state[i][j] ^ key[i + 4 * j]);
				}
			}
			return result;
		}
		inline void AddRoundKeyRuntime(State& state, RoundKeyPart key)
		{
			for (size_t i = 0; i < 4; i++) {
				for (size_t j = 0; j < 4; j++)
				{
					state[i][j] = static_cast<uint8_t>(state[i][j] ^ key[i + 4 * j]);
				}
			}
		}
		inline consteval State SubBytes(State state)
		{
			State result{};
			for (size_t i = 0; i < 4; i++)
			{
				for (size_t j = 0; j < 4; j++)
				{
					uint8_t t = state[i][j];
					result[i][j] = SBox[t / 16][t % 16];
				}
			}
			return result;
		}
		inline void InvSubBytes(State& state) {
			for (size_t i = 0; i < 4; i++)
			{
				for (size_t j = 0; j < 4; j++)
				{
					uint8_t t = state[i][j];
					state[i][j] = invSBox[t / 16][t % 16];
				}
			}
		}
		inline consteval StatePart ShiftRow(StatePart state, int n)    // shift row i on n positions
		{
			StatePart result{};
			for (size_t j = 0; j < 4; j++) {
				result[j] = state[static_cast<size_t>((j + n) % 4)];
			}
			return result;
		}
		inline void ShiftRowRuntime(StatePart& state, int n) {
			unsigned char tmp[4];
			for (int j = 0; j < 4; j++) {
				tmp[j] = state[static_cast<size_t>((j + n) % 4)];
			}
			memcpy(state.data(), tmp, 4 * sizeof(unsigned char));
		}
		inline void InvShiftRows(State& state)
		{
			ShiftRowRuntime(state[1], 4 - 1);
			ShiftRowRuntime(state[2], 4 - 2);
			ShiftRowRuntime(state[3], 4 - 3);
		}
		inline unsigned char mul_bytes(unsigned char a, unsigned char b) // multiplication a and b in galois field
		{
			unsigned char p = 0;
			unsigned char high_bit_mask = 0x80;
			unsigned char high_bit = 0;
			unsigned char modulo = 0x1B; /* x^8 + x^4 + x^3 + x + 1 */


			for (int i = 0; i < 8; i++) {
				if (b & 1) {
					p ^= a;
				}

				high_bit = static_cast<uint8_t>(a & high_bit_mask);
				a <<= 1;
				if (high_bit) {
					a ^= modulo;
				}
				b >>= 1;
			}

			return p;
		}

		inline void InvMixColumns(State& state)
		{
			unsigned char s[4], s1[4];

			for (size_t j = 0; j < 4; j++)
			{
				for (size_t i = 0; i < 4; i++)
				{
					s[i] = state[i][j];
				}
				s1[0] = static_cast<uint8_t>(mul_bytes(0x0e, s[0]) ^ mul_bytes(0x0b, s[1]) ^ mul_bytes(0x0d, s[2]) ^ mul_bytes(0x09, s[3]));
				s1[1] = static_cast<uint8_t>(mul_bytes(0x09, s[0]) ^ mul_bytes(0x0e, s[1]) ^ mul_bytes(0x0b, s[2]) ^ mul_bytes(0x0d, s[3]));
				s1[2] = static_cast<uint8_t>(mul_bytes(0x0d, s[0]) ^ mul_bytes(0x09, s[1]) ^ mul_bytes(0x0e, s[2]) ^ mul_bytes(0x0b, s[3]));
				s1[3] = static_cast<uint8_t>(mul_bytes(0x0b, s[0]) ^ mul_bytes(0x0d, s[1]) ^ mul_bytes(0x09, s[2]) ^ mul_bytes(0x0e, s[3]));

				for (size_t i = 0; i < 4; i++)
				{
					state[i][j] = s1[i];
				}
			}
		}
		inline consteval State ShiftRows(State states) {
			State result{};
			result[0] = states[0];
			result[1] = ShiftRow(states[1], 1);
			result[2] = ShiftRow(states[2], 2);
			result[3] = ShiftRow(states[3], 3);
			return result;
		}
		inline constexpr unsigned char xtime(unsigned char b)
		{
			return static_cast<uint8_t>((b << 1) ^ (((b >> 7) & 1) * 0x1b));
		}
		inline constexpr std::array<uint8_t, 4> MixSingleColumn(std::array<uint8_t, 4> r)
		{
			std::array<uint8_t, 4> result{};
			std::array<uint8_t, 4> a;
			std::array<uint8_t, 4> b;
			uint8_t h;
			for (size_t c = 0; c < 4; c++)
			{
				a[c] = r[c];
				h = (unsigned char)((signed char)r[c] >> 7);
				b[c] = static_cast<uint8_t>(r[c] << 1);
				b[c] ^= 0x1B & h;
			}
			result[0] = static_cast<uint8_t>(b[0] ^ a[3] ^ a[2] ^ b[1] ^ a[1]); /* 2 * a0 + a3 + a2 + 3 * a1 */
			result[1] = static_cast<uint8_t>(b[1] ^ a[0] ^ a[3] ^ b[2] ^ a[2]); /* 2 * a1 + a0 + a3 + 3 * a2 */
			result[2] = static_cast<uint8_t>(b[2] ^ a[1] ^ a[0] ^ b[3] ^ a[3]); /* 2 * a2 + a1 + a0 + 3 * a3 */
			result[3] = static_cast<uint8_t>(b[3] ^ a[2] ^ a[1] ^ b[0] ^ a[0]); /* 2 * a3 + a2 + a1 + 3 * a0 */
			return result;
		}
		inline constexpr State MixColumns(State state)
		{
			State result{};
			std::array<uint8_t, 4> temp;
			for (size_t i = 0; i < 4; ++i)
			{
				for (size_t j = 0; j < 4; ++j)
				{
					temp[j] = state[j][i]; //place the current state column in temp
				}
				temp = MixSingleColumn(temp); //mix it using the wiki implementation
				for (size_t j = 0; j < 4; ++j)
				{
					result[j][i] = temp[j]; //when the column is mixed, place it back into the state
				}
			}
			return result;
		}

		template <size_t keyLength>
		inline consteval Block encryptBlock(Block plain, RoundKey<keyLength> roundKeys) {
			State state{};
			Block result{};
			for (size_t i = 0; i < 4; i++)
				for (size_t j = 0; j < 4; j++)
					state[i][j] = plain[i + j * 4];
			state = AddRoundKey(state, roundKeys[0]);
			for (size_t round = 1; round < parameters<keyLength>::Nr; round++)
			{
				state = SubBytes(state);
				state = ShiftRows(state);
				state = MixColumns(state);
				state = AddRoundKey(state, roundKeys[round]);
			}
			state = SubBytes(state);
			state = ShiftRows(state);
			state = AddRoundKey(state, roundKeys[parameters<keyLength>::Nr]);
			for (size_t i = 0; i < 4; i++)
				for (size_t j = 0; j < 4; j++)
					result[i + j * 4] = state[i][j];
			return result;
		}
		template <size_t keyLength>
		inline void DecryptBlock(Block crypted, RoundKey<keyLength> roundKeys, uint8_t* result)
		{
			State state{};

			for (size_t i = 0; i < 4; i++)
			{
				for (size_t j = 0; j < 4; j++) {
					state[i][j] = crypted[i + 4 * j];
				}
			}

			AddRoundKeyRuntime(state, roundKeys[parameters<keyLength>::Nr]);

			for (size_t round = parameters<keyLength>::Nr - 1; round >= 1; round--)
			{
				InvSubBytes(state);
				InvShiftRows(state);
				AddRoundKeyRuntime(state, roundKeys[round]);
				InvMixColumns(state);
			}
			InvSubBytes(state);
			InvShiftRows(state);
			AddRoundKeyRuntime(state, roundKeys[0]);

			for (size_t i = 0; i < 4; i++)
			{
				for (size_t j = 0; j < 4; j++) {
					result[i + 4 * j] = state[i][j];
				}
			}
		}

		template <size_t keyLength>
		using Key = std::array<uint8_t, 4 * parameters<keyLength>::Nk>;


		inline consteval std::array<uint8_t, 4> RotWord(std::array<uint8_t, 4> x)
		{
			return std::array<uint8_t, 4>{x[1], x[2], x[3], x[0]};
		}

		inline consteval std::array<uint8_t, 4> SubWord(std::array<uint8_t, 4> x)
		{
			std::array<uint8_t, 4> result{};
			for (size_t i = 0; i < 4; i++)
			{
				result[i] = SBox[x[i] / 16][x[i] % 16];
			}
			return result;
		}
		inline void RotWordInplace(std::array<uint8_t, 4>& x)
		{
			uint8_t tmp = x[0];
			x[0] = x[1];
			x[1] = x[2];
			x[2] = x[3];
			x[3] = tmp;
		}
		inline void SubWordInplace(std::array<uint8_t, 4>& x)
		{
			for (size_t i = 0; i < 4; i++)
				x[i] = SBox[x[i] / 16][x[i] % 16];
		}
		inline consteval std::array<uint8_t, 4> XorWords(std::array<uint8_t, 4>a, std::array<uint8_t, 4>b)
		{
			std::array<uint8_t, 4> result{};
			for (size_t i = 0; i < 4; i++)
			{
				result[i] = static_cast<uint8_t>(a[i] ^ b[i]);
			}
			return result;
		}
		inline void XorWordsInplace(std::array<uint8_t, 4>& a, std::array<uint8_t, 4>b)
		{
			for (size_t i = 0; i < 4; i++)
				a[i] ^= b[i];
		}
		inline constexpr std::array<uint8_t, 4> Rcon(size_t n)
		{
			std::array<uint8_t, 4> result{};
			unsigned char c = 1;
			for (size_t i = 0; i < n - 1; i++)
			{
				c = xtime(c);
			}

			result[0] = c;
			result[1] = result[2] = result[3] = 0;
			return result;
		}

		template <size_t keyLength>
		inline consteval RoundKey<keyLength> keyExpansion(Key<keyLength> key)
		{
			RoundKey<keyLength> result{};
			std::array<uint8_t, 4> temp;
			for (size_t i = 0; i < key.size(); i++)
				result[i / BlockSize][i % BlockSize] = key[i];
			for (size_t i = key.size(); i < BlockSize * result.size(); i += 4) {
				temp[0] = result[(i - 4 + 0) / BlockSize][(i - 4 + 0) % BlockSize];
				temp[1] = result[(i - 4 + 1) / BlockSize][(i - 4 + 1) % BlockSize];
				temp[2] = result[(i - 4 + 2) / BlockSize][(i - 4 + 2) % BlockSize];
				temp[3] = result[(i - 4 + 3) / BlockSize][(i - 4 + 3) % BlockSize];
				if (i / 4 % parameters<keyLength>::Nk == 0)
				{
					temp = XorWords(SubWord(RotWord(temp)), Rcon(i / (parameters<keyLength>::Nk * 4)));
				}
				else if (parameters<keyLength>::Nk > 6 && i / 4 % parameters<keyLength>::Nk == 4)
				{
					temp = SubWord(temp);
				}
				result[(i + 0) / BlockSize][(i + 0) % BlockSize] =
					static_cast<uint8_t>(
						result[(i + 0 - 4 * parameters<keyLength>::Nk) / BlockSize][(i + 0 - 4 * parameters<keyLength>::Nk) % BlockSize] ^ temp[0]
						);
				result[(i + 1) / BlockSize][(i + 1) % BlockSize] =
					static_cast<uint8_t>(
						result[(i + 1 - 4 * parameters<keyLength>::Nk) / BlockSize][(i + 1 - 4 * parameters<keyLength>::Nk) % BlockSize] ^ temp[1]
						);
				result[(i + 2) / BlockSize][(i + 2) % BlockSize] =
					static_cast<uint8_t>(
						result[(i + 2 - 4 * parameters<keyLength>::Nk) / BlockSize][(i + 2 - 4 * parameters<keyLength>::Nk) % BlockSize] ^ temp[2]
						);
				result[(i + 3) / BlockSize][(i + 3) % BlockSize] =
					static_cast<uint8_t>(
						result[(i + 3 - 4 * parameters<keyLength>::Nk) / BlockSize][(i + 3 - 4 * parameters<keyLength>::Nk) % BlockSize] ^ temp[3]
						);
			}
			return result;
		}
		template <size_t keyLength>
		inline void keyExpansionRuntime(Key<keyLength> key, RoundKey<keyLength>& roundKey)
		{
			std::array<uint8_t, 4> temp;
			for (size_t i = 0; i < key.size(); i++)
				roundKey[i / BlockSize][i % BlockSize] = key[i];
			for (size_t i = key.size(); i < BlockSize * roundKey.size(); i += 4) {
				temp[0] = roundKey[(i - 4 + 0) / BlockSize][(i - 4 + 0) % BlockSize];
				temp[1] = roundKey[(i - 4 + 1) / BlockSize][(i - 4 + 1) % BlockSize];
				temp[2] = roundKey[(i - 4 + 2) / BlockSize][(i - 4 + 2) % BlockSize];
				temp[3] = roundKey[(i - 4 + 3) / BlockSize][(i - 4 + 3) % BlockSize];
				if (i / 4 % parameters<keyLength>::Nk == 0)
				{
					RotWordInplace(temp);
					SubWordInplace(temp);
					XorWordsInplace(temp, Rcon(i / (parameters<keyLength>::Nk * 4)));
				}
				else if (parameters<keyLength>::Nk > 6 && i / 4 % parameters<keyLength>::Nk == 4)
				{
					SubWordInplace(temp);
				}
				roundKey[(i + 0) / BlockSize][(i + 0) % BlockSize] =
					static_cast<uint8_t>(roundKey[(i + 0 - 4 * parameters<keyLength>::Nk) / BlockSize][(i + 0 - 4 * parameters<keyLength>::Nk) % BlockSize] ^ temp[0]);
				roundKey[(i + 1) / BlockSize][(i + 1) % BlockSize] =
					static_cast<uint8_t>(roundKey[(i + 1 - 4 * parameters<keyLength>::Nk) / BlockSize][(i + 1 - 4 * parameters<keyLength>::Nk) % BlockSize] ^ temp[1]);
				roundKey[(i + 2) / BlockSize][(i + 2) % BlockSize] =
					static_cast<uint8_t>(roundKey[(i + 2 - 4 * parameters<keyLength>::Nk) / BlockSize][(i + 2 - 4 * parameters<keyLength>::Nk) % BlockSize] ^ temp[2]);
				roundKey[(i + 3) / BlockSize][(i + 3) % BlockSize] =
					static_cast<uint8_t>(roundKey[(i + 3 - 4 * parameters<keyLength>::Nk) / BlockSize][(i + 3 - 4 * parameters<keyLength>::Nk) % BlockSize] ^ temp[3]);
			}
		}
		template <size_t inputLen>
		inline consteval std::array<std::array<uint8_t, BlockSize>, (inputLen - 1) / BlockSize + 1>
			alignInput(std::array<uint8_t, inputLen> input)
		{
			std::array<std::array<uint8_t, BlockSize>, (inputLen - 1) / BlockSize + 1> result{};
			for (size_t i = 0; i < inputLen / BlockSize; i++)
				for (size_t j = 0; j < BlockSize; j++)
					result[i][j] = input[i * BlockSize + j];
			if constexpr (inputLen % BlockSize != 0) {
				for (size_t j = 0; j < inputLen % BlockSize; j++) {
					result[inputLen / BlockSize][j] = input[inputLen - inputLen % BlockSize + j];
				}
				for (size_t j = inputLen % BlockSize; j < BlockSize; j++)
					result[inputLen / BlockSize][j] = 0;
			}
			return result;
		}
		template <size_t keyLength, size_t inputLen>
		inline consteval auto encryptCBC(std::array<uint8_t, inputLen> input, Key<keyLength> key, Block iv) {
			auto aligned_input = alignInput(input);
			decltype(aligned_input) result{};
			for (size_t i = 0; i < aligned_input.size(); i++) {
				result[i] = aligned_input[i];
			}
			auto roundKeys = keyExpansion<keyLength>(key);
			auto block = iv;
			for (size_t i = 0; i < aligned_input.size(); i++)
			{
				block = xorBlock(aligned_input[i], block);
				result[i] = encryptBlock<256>(block, roundKeys);
				block = result[i];
			}
			return result;
		}
		template <size_t keyLength, size_t inputLen>
		inline void DecryptCBC(
			std::array<std::array<uint8_t, BlockSize>, (inputLen - 1) / BlockSize + 1> encrypted,
			std::array<uint8_t, inputLen>& output,
			Block iv,
			Key<keyLength> key
		)
		{
			RoundKey<keyLength> roundKey;
			keyExpansionRuntime<keyLength>(key, roundKey);
			auto block = iv;
			for (size_t i = 0; i < encrypted.size(); i++)
			{
				DecryptBlock<keyLength>(encrypted[i], roundKey, output.data() + i * BlockSize);
				xorBlockInplace(output.data() + i * BlockSize, block);
				block = encrypted[i];
			}
		}
		template<uint64_t id>
		class randomBlock {
			inline static constexpr uint64_t id2 = id + 0x0101010101010101ull;
			inline static constexpr std::array<uint8_t, 16> plain = {
				static_cast<uint8_t>(id),
				static_cast<uint8_t>(id / 256),
				static_cast<uint8_t>(id / 256 / 256),
				static_cast<uint8_t>(id / 256 / 256 / 256),
				static_cast<uint8_t>(id / 256 / 256 / 256 / 256),
				static_cast<uint8_t>(id / 256 / 256 / 256 / 256 / 256),
				static_cast<uint8_t>(id / 256 / 256 / 256 / 256 / 256 / 256),
				static_cast<uint8_t>(id / 256 / 256 / 256 / 256 / 256 / 256 / 256),
				static_cast<uint8_t>(id2),
				static_cast<uint8_t>(id2 / 256),
				static_cast<uint8_t>(id2 / 256 / 256),
				static_cast<uint8_t>(id2 / 256 / 256 / 256),
				static_cast<uint8_t>(id2 / 256 / 256 / 256 / 256),
				static_cast<uint8_t>(id2 / 256 / 256 / 256 / 256 / 256),
				static_cast<uint8_t>(id2 / 256 / 256 / 256 / 256 / 256 / 256),
				static_cast<uint8_t>(id2 / 256 / 256 / 256 / 256 / 256 / 256 / 256),
			};

			inline static constexpr uint64_t id3 = id + 0x0110011001100110ull;
			inline static constexpr uint64_t id4 = id + 0x1100010111000101ull;
			inline static constexpr std::array<uint8_t, 32> rKey = {
				static_cast<uint8_t>(id3),
				static_cast<uint8_t>(id3 / 256),
				static_cast<uint8_t>(id3 / 256 / 256),
				static_cast<uint8_t>(id3 / 256 / 256 / 256),
				static_cast<uint8_t>(id3 / 256 / 256 / 256 / 256),
				static_cast<uint8_t>(id3 / 256 / 256 / 256 / 256 / 256),
				static_cast<uint8_t>(id3 / 256 / 256 / 256 / 256 / 256 / 256),
				static_cast<uint8_t>(id3 / 256 / 256 / 256 / 256 / 256 / 256 / 256),
				static_cast<uint8_t>(id3 / 256 / 256 / 256 / 256 / 256 / 256 / 256),
				static_cast<uint8_t>(id3 / 256 / 256 / 256 / 256 / 256 / 256),
				static_cast<uint8_t>(id3 / 256 / 256 / 256 / 256 / 256),
				static_cast<uint8_t>(id3 / 256 / 256 / 256 / 256),
				static_cast<uint8_t>(id3 / 256 / 256 / 256),
				static_cast<uint8_t>(id3 / 256 / 256),
				static_cast<uint8_t>(id3 / 256),
				static_cast<uint8_t>(id3),
				static_cast<uint8_t>(id4),
				static_cast<uint8_t>(id4 / 256),
				static_cast<uint8_t>(id4 / 256 / 256),
				static_cast<uint8_t>(id4 / 256 / 256 / 256),
				static_cast<uint8_t>(id4 / 256 / 256 / 256 / 256),
				static_cast<uint8_t>(id4 / 256 / 256 / 256 / 256 / 256),
				static_cast<uint8_t>(id4 / 256 / 256 / 256 / 256 / 256 / 256),
				static_cast<uint8_t>(id4 / 256 / 256 / 256 / 256 / 256 / 256 / 256),
				static_cast<uint8_t>(id4 / 256 / 256 / 256 / 256 / 256 / 256 / 256),
				static_cast<uint8_t>(id4 / 256 / 256 / 256 / 256 / 256 / 256),
				static_cast<uint8_t>(id4 / 256 / 256 / 256 / 256 / 256),
				static_cast<uint8_t>(id4 / 256 / 256 / 256 / 256),
				static_cast<uint8_t>(id4 / 256 / 256 / 256),
				static_cast<uint8_t>(id4 / 256 / 256),
				static_cast<uint8_t>(id4 / 256),
				static_cast<uint8_t>(id4)
			};
		public:
			constexpr static std::array<uint8_t, 16> data = encryptCBC<256, 16>(plain, rKey, plain)[0];
		};

		inline consteval uint64_t u64_from_block(std::array<uint8_t, 16> x) {
			uint64_t tmp = 0;
			for (size_t i = 0; i < 8; i++) {
				tmp |= static_cast<uint64_t>(x[i] ^ x[(i * 3) % 8 + 8]) << (i * 8);
			}
			return tmp;
		}
		template<uint64_t id, size_t N>
		inline consteval std::array<uint8_t, 16 * N> genRandomBlocks() {
			std::array<uint8_t, 16 * N> result{};
			constexpr std::array<uint8_t, 16> block1 = randomBlock<id>::data;
			if constexpr (N == 1) {
				result = block1;
			}
			else {
				constexpr std::array<uint8_t, 16> block2 = randomBlock<id ^ 0x1111111100000000ull>::data;
				constexpr uint64_t next_id = u64_from_block(block2);
				constexpr auto left = genRandomBlocks<next_id, N - 1>();
				for (size_t i = 0; i < 16 * (N - 1); i++) {
					result[i] = left[i];
				}
				for (size_t i = 0; i < 16; i++) {
					result[i + 16 * (N - 1)] = block1[i];
				}
			}
			return result;
		}
		template<uint64_t id, size_t N>
		class randomBlocks {
		public:
			constexpr static std::array<uint8_t, 16 * N> data = genRandomBlocks<id, N>();
		};
		template <typename T, size_t N>
		class ProtectedMemory {
		public:
			inline ProtectedMemory() {};
			std::array<uint8_t, sizeof(T)* N> plain;
			OBFUSCATOR_INLINE const T* data() {
				return reinterpret_cast<const T*>(plain.data());
			}
			OBFUSCATOR_INLINE void clear() {
				memset(plain.data(), 0, plain.size());
			}
			inline ~ProtectedMemory() {
				clear();
			}
		};

		template<uint64_t seed, size_t N>
		class ProtectedResource {
		private:
			ProtectedResource& operator=(const ProtectedResource& other) = delete;
		public:
			inline consteval ProtectedResource(const ProtectedResource<seed, N>& other) :encrypted(other.encrypted)
			{}
			const std::array<std::array<uint8_t, BlockSize>, (N - 1) / BlockSize + 1> encrypted;
			constexpr static std::array<uint8_t, 32> key = randomBlocks<seed, 2>::data;
			constexpr static std::array<uint8_t, 16> iv = randomBlocks<
				u64_from_block(randomBlock<seed ^ 0x0011011011001001ull>::data)
				, 1>::data;
			inline consteval ProtectedResource(std::array<uint8_t, N> data) :
				encrypted(encryptCBC<256, N>(data, key, iv))
			{}
			enum { size = N };
		};
		template<uint64_t seed, typename T, size_t N>
		OBFUSCATOR_INLINE ProtectedMemory<T, N> decryptResource(const ProtectedResource<seed, sizeof(T)* N>& x) {
			ProtectedMemory<T, N> result;
			DecryptCBC<256, sizeof(T)* N>(x.encrypted, result.plain, x.iv, x.key);
			return result;
		}
		template <typename T, size_t N>
		inline consteval std::array<uint8_t, sizeof(T)* N> convert(std::array<T, N > x) {
			struct {
				uint8_t data[sizeof(T) * N];
			} target;
			struct {
				T data[N];
			} source;
			for (size_t i = 0; i < N; i++)
				source.data[i] = x[i];
			target = std::bit_cast<decltype(target)>(source);
			std::array<uint8_t, sizeof(T)* N> result{};
			for (size_t i = 0; i < sizeof(T) * N; i++)
				result[i] = target.data[i];
			return result;
		};
		template <typename T, size_t N>
		inline consteval std::array<T, N> make_array(const T data[N]) {
			std::array<T, N> result{};
			for (size_t i = 0; i < N; i++)
				result[i] = data[i];
			return result;
		}
		template <uint64_t seed, typename T, size_t N>
		class ProtectedArray :public ProtectedResource<seed, sizeof(T)* N> {
		private:
			ProtectedArray& operator=(const ProtectedArray<seed, T, N>& other) = delete;
		public:
			inline consteval ProtectedArray(const ProtectedArray<seed, T, N>& other) = delete;
			inline consteval ProtectedArray(std::array<T, N> data) :
				ProtectedResource<seed, sizeof(T)* N>(convert(data)) {}
			inline consteval ProtectedArray(const T data[N]) : ProtectedArray(make_array<T, N>(data))
			{}
			inline ProtectedMemory<T, N> decrypt()const {
				return decryptResource<seed, T, N>(*this);
			}
		};

		template<size_t filenameSize>
		consteval uint64_t genRandom(uint64_t line, const char* filename, const char* date, const char* time) {
			std::array<uint8_t, sizeof(uint64_t) + 8 + 11 + filenameSize> input;
			for (size_t i = 0; i < sizeof(uint64_t); i++)
				input[i] = static_cast<uint8_t>(line >> (i * 8));
			for (size_t i = 0; i < filenameSize; i++)
				input[i + sizeof(uint64_t)] = static_cast<uint8_t>(filename[i]);
			for (size_t i = 0; i < 11; i++)
				input[i + sizeof(uint64_t) + filenameSize] = static_cast<uint8_t>(date[i]);
			for (size_t i = 0; i < 8; i++)
				input[i + sizeof(uint64_t) + filenameSize + 11] = static_cast<uint8_t>(time[i]);
			Block iv;
			Key<256> key;
			for (size_t i = 0; i < iv.size(); i++)
				iv[i] = 0;
			for (size_t i = 0; i < key.size(); i++)
				key[i] = 0;
			for (size_t i = 0; i < input.size(); i++) {
				iv[i % iv.size()] ^= input[i];
				key[i % key.size()] ^= input[i];
			}
			return u64_from_block(encryptCBC<256, input.size()>(input, key, iv)[0]);
		}
		consteval int c_strlen(const char* str)
		{
			return *str ? 1 + c_strlen(str + 1) : 0;
		}
	}
	using Impl::ProtectedMemory;
	using Impl::ProtectedResource;
	using Impl::ProtectedArray;
#undef OBFUSCATOR_INLINE
#ifndef OBFUSCATOR_NO_MACROS
#define RANDOM_SEED	(SimpleObfuscator::Impl::genRandom<SimpleObfuscator::Impl::c_strlen(__FILE__)>(__LINE__,__FILE__,__DATE__,__TIME__))
#define DECLARE_PROTECTED(type,count,name,...) constexpr SimpleObfuscator::ProtectedArray<RANDOM_SEED,type,count> name(__VA_ARGS__);
#define PROTECTED(type,count,...) SimpleObfuscator::ProtectedArray<RANDOM_SEED,type,count>(__VA_ARGS__).decrypt().data()
#define PROTECTED_STRING(...) PROTECTED(char,SimpleObfuscator::Impl::c_strlen(__VA_ARGS__),__VA_ARGS__)
#endif
};
