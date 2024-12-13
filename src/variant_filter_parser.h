// A Bison parser, made by GNU Bison 3.8.2.

// Skeleton interface for Bison LALR(1) parsers in C++

// Copyright (C) 2002-2015, 2018-2021 Free Software Foundation, Inc.

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

// As a special exception, you may create a larger work that contains
// part or all of the Bison parser skeleton and distribute that work
// under terms of your choice, so long as that work isn't itself a
// parser generator using the skeleton or a modified version thereof
// as a parser skeleton.  Alternatively, if you modify or redistribute
// the parser skeleton itself, you may (at your option) remove this
// special exception, which will cause the skeleton and the resulting
// Bison output files to be licensed under the GNU General Public
// License without this special exception.

// This special exception was added by the Free Software Foundation in
// version 2.2 of Bison.


/**
 ** \file variant_filter_parser.h
 ** Define the kim::parser class.
 */

// C++ LALR(1) parser skeleton written by Akim Demaille.

// DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
// especially those whose name start with YY_ or yy_.  They are
// private implementation details that can be changed or removed.

#ifndef YY_YY_VARIANT_FILTER_PARSER_H_INCLUDED
# define YY_YY_VARIANT_FILTER_PARSER_H_INCLUDED
// "%code requires" blocks.
#line 21 "variant_filter_parser.ypp"

#include <variant_filter.h>
namespace kim {
  class VariantFilterDriver;
  /**
   * \class VariantFilterParser::context
   *
   *  The Context of VariantFilterParser.
   */
}

#line 61 "variant_filter_parser.h"

# include <cassert>
# include <cstdlib> // std::abort
# include <iostream>
# include <stdexcept>
# include <string>
# include <vector>

#if defined __cplusplus
# define YY_CPLUSPLUS __cplusplus
#else
# define YY_CPLUSPLUS 199711L
#endif

// Support move semantics when possible.
#if 201103L <= YY_CPLUSPLUS
# define YY_MOVE           std::move
# define YY_MOVE_OR_COPY   move
# define YY_MOVE_REF(Type) Type&&
# define YY_RVREF(Type)    Type&&
# define YY_COPY(Type)     Type
#else
# define YY_MOVE
# define YY_MOVE_OR_COPY   copy
# define YY_MOVE_REF(Type) Type&
# define YY_RVREF(Type)    const Type&
# define YY_COPY(Type)     const Type&
#endif

// Support noexcept when possible.
#if 201103L <= YY_CPLUSPLUS
# define YY_NOEXCEPT noexcept
# define YY_NOTHROW
#else
# define YY_NOEXCEPT
# define YY_NOTHROW throw ()
#endif

// Support constexpr when possible.
#if 201703 <= YY_CPLUSPLUS
# define YY_CONSTEXPR constexpr
#else
# define YY_CONSTEXPR
#endif

#include <typeinfo>
#ifndef YY_ASSERT
# include <cassert>
# define YY_ASSERT assert
#endif


#ifndef YY_ATTRIBUTE_PURE
# if defined __GNUC__ && 2 < __GNUC__ + (96 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_PURE __attribute__ ((__pure__))
# else
#  define YY_ATTRIBUTE_PURE
# endif
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# if defined __GNUC__ && 2 < __GNUC__ + (7 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_UNUSED __attribute__ ((__unused__))
# else
#  define YY_ATTRIBUTE_UNUSED
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YY_USE(E) ((void) (E))
#else
# define YY_USE(E) /* empty */
#endif

/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
#if defined __GNUC__ && ! defined __ICC && 406 <= __GNUC__ * 100 + __GNUC_MINOR__
# if __GNUC__ * 100 + __GNUC_MINOR__ < 407
#  define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN                           \
    _Pragma ("GCC diagnostic push")                                     \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")
# else
#  define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN                           \
    _Pragma ("GCC diagnostic push")                                     \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")              \
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# endif
# define YY_IGNORE_MAYBE_UNINITIALIZED_END      \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif

#if defined __cplusplus && defined __GNUC__ && ! defined __ICC && 6 <= __GNUC__
# define YY_IGNORE_USELESS_CAST_BEGIN                          \
    _Pragma ("GCC diagnostic push")                            \
    _Pragma ("GCC diagnostic ignored \"-Wuseless-cast\"")
# define YY_IGNORE_USELESS_CAST_END            \
    _Pragma ("GCC diagnostic pop")
#endif
#ifndef YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_END
#endif

# ifndef YY_CAST
#  ifdef __cplusplus
#   define YY_CAST(Type, Val) static_cast<Type> (Val)
#   define YY_REINTERPRET_CAST(Type, Val) reinterpret_cast<Type> (Val)
#  else
#   define YY_CAST(Type, Val) ((Type) (Val))
#   define YY_REINTERPRET_CAST(Type, Val) ((Type) (Val))
#  endif
# endif
# ifndef YY_NULLPTR
#  if defined __cplusplus
#   if 201103L <= __cplusplus
#    define YY_NULLPTR nullptr
#   else
#    define YY_NULLPTR 0
#   endif
#  else
#   define YY_NULLPTR ((void*)0)
#  endif
# endif

/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 1
#endif

#line 15 "variant_filter_parser.ypp"
namespace kim {
#line 202 "variant_filter_parser.h"


  /// A point in a source file.
  class position
  {
  public:
    /// Type for file name.
    typedef const std::string filename_type;
    /// Type for line and column numbers.
    typedef int counter_type;

    /// Construct a position.
    explicit position (filename_type* f = YY_NULLPTR,
                       counter_type l = 1,
                       counter_type c = 1)
      : filename (f)
      , line (l)
      , column (c)
    {}


    /// Initialization.
    void initialize (filename_type* fn = YY_NULLPTR,
                     counter_type l = 1,
                     counter_type c = 1)
    {
      filename = fn;
      line = l;
      column = c;
    }

    /** \name Line and Column related manipulators
     ** \{ */
    /// (line related) Advance to the COUNT next lines.
    void lines (counter_type count = 1)
    {
      if (count)
        {
          column = 1;
          line = add_ (line, count, 1);
        }
    }

    /// (column related) Advance to the COUNT next columns.
    void columns (counter_type count = 1)
    {
      column = add_ (column, count, 1);
    }
    /** \} */

    /// File name to which this position refers.
    filename_type* filename;
    /// Current line number.
    counter_type line;
    /// Current column number.
    counter_type column;

  private:
    /// Compute max (min, lhs+rhs).
    static counter_type add_ (counter_type lhs, counter_type rhs, counter_type min)
    {
      return lhs + rhs < min ? min : lhs + rhs;
    }
  };

  /// Add \a width columns, in place.
  inline position&
  operator+= (position& res, position::counter_type width)
  {
    res.columns (width);
    return res;
  }

  /// Add \a width columns.
  inline position
  operator+ (position res, position::counter_type width)
  {
    return res += width;
  }

  /// Subtract \a width columns, in place.
  inline position&
  operator-= (position& res, position::counter_type width)
  {
    return res += -width;
  }

  /// Subtract \a width columns.
  inline position
  operator- (position res, position::counter_type width)
  {
    return res -= width;
  }

  /** \brief Intercept output stream redirection.
   ** \param ostr the destination output stream
   ** \param pos a reference to the position to redirect
   */
  template <typename YYChar>
  std::basic_ostream<YYChar>&
  operator<< (std::basic_ostream<YYChar>& ostr, const position& pos)
  {
    if (pos.filename)
      ostr << *pos.filename << ':';
    return ostr << pos.line << '.' << pos.column;
  }

  /// Two points in a source file.
  class location
  {
  public:
    /// Type for file name.
    typedef position::filename_type filename_type;
    /// Type for line and column numbers.
    typedef position::counter_type counter_type;

    /// Construct a location from \a b to \a e.
    location (const position& b, const position& e)
      : begin (b)
      , end (e)
    {}

    /// Construct a 0-width location in \a p.
    explicit location (const position& p = position ())
      : begin (p)
      , end (p)
    {}

    /// Construct a 0-width location in \a f, \a l, \a c.
    explicit location (filename_type* f,
                       counter_type l = 1,
                       counter_type c = 1)
      : begin (f, l, c)
      , end (f, l, c)
    {}


    /// Initialization.
    void initialize (filename_type* f = YY_NULLPTR,
                     counter_type l = 1,
                     counter_type c = 1)
    {
      begin.initialize (f, l, c);
      end = begin;
    }

    /** \name Line and Column related manipulators
     ** \{ */
  public:
    /// Reset initial location to final location.
    void step ()
    {
      begin = end;
    }

    /// Extend the current location to the COUNT next columns.
    void columns (counter_type count = 1)
    {
      end += count;
    }

    /// Extend the current location to the COUNT next lines.
    void lines (counter_type count = 1)
    {
      end.lines (count);
    }
    /** \} */


  public:
    /// Beginning of the located region.
    position begin;
    /// End of the located region.
    position end;
  };

  /// Join two locations, in place.
  inline location&
  operator+= (location& res, const location& end)
  {
    res.end = end.end;
    return res;
  }

  /// Join two locations.
  inline location
  operator+ (location res, const location& end)
  {
    return res += end;
  }

  /// Add \a width columns to the end position, in place.
  inline location&
  operator+= (location& res, location::counter_type width)
  {
    res.columns (width);
    return res;
  }

  /// Add \a width columns to the end position.
  inline location
  operator+ (location res, location::counter_type width)
  {
    return res += width;
  }

  /// Subtract \a width columns to the end position, in place.
  inline location&
  operator-= (location& res, location::counter_type width)
  {
    return res += -width;
  }

  /// Subtract \a width columns to the end position.
  inline location
  operator- (location res, location::counter_type width)
  {
    return res -= width;
  }

  /** \brief Intercept output stream redirection.
   ** \param ostr the destination output stream
   ** \param loc a reference to the location to redirect
   **
   ** Avoid duplicate information.
   */
  template <typename YYChar>
  std::basic_ostream<YYChar>&
  operator<< (std::basic_ostream<YYChar>& ostr, const location& loc)
  {
    location::counter_type end_col
      = 0 < loc.end.column ? loc.end.column - 1 : 0;
    ostr << loc.begin;
    if (loc.end.filename
        && (!loc.begin.filename
            || *loc.begin.filename != *loc.end.filename))
      ostr << '-' << loc.end.filename << ':' << loc.end.line << '.' << end_col;
    else if (loc.begin.line < loc.end.line)
      ostr << '-' << loc.end.line << '.' << end_col;
    else if (loc.begin.column < end_col)
      ostr << '-' << end_col;
    return ostr;
  }


  /// A Bison parser.
  class VariantFilterParser
  {
  public:
#ifdef YYSTYPE
# ifdef __GNUC__
#  pragma GCC message "bison: do not #define YYSTYPE in C++, use %define api.value.type"
# endif
    typedef YYSTYPE value_type;
#else
  /// A buffer to store and retrieve objects.
  ///
  /// Sort of a variant, but does not keep track of the nature
  /// of the stored data, since that knowledge is available
  /// via the current parser state.
  class value_type
  {
  public:
    /// Type of *this.
    typedef value_type self_type;

    /// Empty construction.
    value_type () YY_NOEXCEPT
      : yyraw_ ()
      , yytypeid_ (YY_NULLPTR)
    {}

    /// Construct and fill.
    template <typename T>
    value_type (YY_RVREF (T) t)
      : yytypeid_ (&typeid (T))
    {
      YY_ASSERT (sizeof (T) <= size);
      new (yyas_<T> ()) T (YY_MOVE (t));
    }

#if 201103L <= YY_CPLUSPLUS
    /// Non copyable.
    value_type (const self_type&) = delete;
    /// Non copyable.
    self_type& operator= (const self_type&) = delete;
#endif

    /// Destruction, allowed only if empty.
    ~value_type () YY_NOEXCEPT
    {
      YY_ASSERT (!yytypeid_);
    }

# if 201103L <= YY_CPLUSPLUS
    /// Instantiate a \a T in here from \a t.
    template <typename T, typename... U>
    T&
    emplace (U&&... u)
    {
      YY_ASSERT (!yytypeid_);
      YY_ASSERT (sizeof (T) <= size);
      yytypeid_ = & typeid (T);
      return *new (yyas_<T> ()) T (std::forward <U>(u)...);
    }
# else
    /// Instantiate an empty \a T in here.
    template <typename T>
    T&
    emplace ()
    {
      YY_ASSERT (!yytypeid_);
      YY_ASSERT (sizeof (T) <= size);
      yytypeid_ = & typeid (T);
      return *new (yyas_<T> ()) T ();
    }

    /// Instantiate a \a T in here from \a t.
    template <typename T>
    T&
    emplace (const T& t)
    {
      YY_ASSERT (!yytypeid_);
      YY_ASSERT (sizeof (T) <= size);
      yytypeid_ = & typeid (T);
      return *new (yyas_<T> ()) T (t);
    }
# endif

    /// Instantiate an empty \a T in here.
    /// Obsolete, use emplace.
    template <typename T>
    T&
    build ()
    {
      return emplace<T> ();
    }

    /// Instantiate a \a T in here from \a t.
    /// Obsolete, use emplace.
    template <typename T>
    T&
    build (const T& t)
    {
      return emplace<T> (t);
    }

    /// Accessor to a built \a T.
    template <typename T>
    T&
    as () YY_NOEXCEPT
    {
      YY_ASSERT (yytypeid_);
      YY_ASSERT (*yytypeid_ == typeid (T));
      YY_ASSERT (sizeof (T) <= size);
      return *yyas_<T> ();
    }

    /// Const accessor to a built \a T (for %printer).
    template <typename T>
    const T&
    as () const YY_NOEXCEPT
    {
      YY_ASSERT (yytypeid_);
      YY_ASSERT (*yytypeid_ == typeid (T));
      YY_ASSERT (sizeof (T) <= size);
      return *yyas_<T> ();
    }

    /// Swap the content with \a that, of same type.
    ///
    /// Both variants must be built beforehand, because swapping the actual
    /// data requires reading it (with as()), and this is not possible on
    /// unconstructed variants: it would require some dynamic testing, which
    /// should not be the variant's responsibility.
    /// Swapping between built and (possibly) non-built is done with
    /// self_type::move ().
    template <typename T>
    void
    swap (self_type& that) YY_NOEXCEPT
    {
      YY_ASSERT (yytypeid_);
      YY_ASSERT (*yytypeid_ == *that.yytypeid_);
      std::swap (as<T> (), that.as<T> ());
    }

    /// Move the content of \a that to this.
    ///
    /// Destroys \a that.
    template <typename T>
    void
    move (self_type& that)
    {
# if 201103L <= YY_CPLUSPLUS
      emplace<T> (std::move (that.as<T> ()));
# else
      emplace<T> ();
      swap<T> (that);
# endif
      that.destroy<T> ();
    }

# if 201103L <= YY_CPLUSPLUS
    /// Move the content of \a that to this.
    template <typename T>
    void
    move (self_type&& that)
    {
      emplace<T> (std::move (that.as<T> ()));
      that.destroy<T> ();
    }
#endif

    /// Copy the content of \a that to this.
    template <typename T>
    void
    copy (const self_type& that)
    {
      emplace<T> (that.as<T> ());
    }

    /// Destroy the stored \a T.
    template <typename T>
    void
    destroy ()
    {
      as<T> ().~T ();
      yytypeid_ = YY_NULLPTR;
    }

  private:
#if YY_CPLUSPLUS < 201103L
    /// Non copyable.
    value_type (const self_type&);
    /// Non copyable.
    self_type& operator= (const self_type&);
#endif

    /// Accessor to raw memory as \a T.
    template <typename T>
    T*
    yyas_ () YY_NOEXCEPT
    {
      void *yyp = yyraw_;
      return static_cast<T*> (yyp);
     }

    /// Const accessor to raw memory as \a T.
    template <typename T>
    const T*
    yyas_ () const YY_NOEXCEPT
    {
      const void *yyp = yyraw_;
      return static_cast<const T*> (yyp);
     }

    /// An auxiliary type to compute the largest semantic type.
    union union_type
    {
      // "boolean value"
      // filter
      // simple_expr
      // boolean
      char dummy1[sizeof (bool)];

      // "number value"
      // number
      char dummy2[sizeof (double)];

      // "integer value"
      // integer
      char dummy3[sizeof (int)];

      // bool_op
      char dummy4[sizeof (kim::VariantFilter::BooleanOperator)];

      // "equality operator ('=')"
      // "regular expression match operator ('~')"
      // "inequality operator ('!=' or '<>')"
      // "order relationship operator '<'"
      // "order relationship operator '<='"
      // "order relationship operator '>'"
      // "order relationship operator '>='"
      // any_op
      char dummy5[sizeof (kim::VariantFilter::ComparisonOperator)];

      // "conjunction operator ('AND' or '&&')"
      // "disjunction operator ('or' or '||')"
      // "negation operator ('not' or '!')"
      // "opening parenthesis '('"
      // "closing parenthesis ')'"
      char dummy6[sizeof (kim::VariantFilter::ExpressionOperator)];

      // "field attribute 'CHROM'"
      // "field attribute 'POS'"
      // "field attribute 'ID'"
      // "field attribute 'QUAL'"
      // "field attribute 'FILTER'"
      // "field attribute 'INFO'"
      // "status 'Ploidy'"
      // "status 'SNP' (Single Nucleotide Polymorphism)"
      // "status 'MSNP' (multi-allelic SNP)"
      // "status 'SV' (Structural Variant)"
      // "status 'InDel' (Insertion-Deletion)"
      char dummy7[sizeof (kim::VariantFilter::FilterOperation)];

      // num_op
      char dummy8[sizeof (kim::VariantFilter::NumericalOperator)];

      // str_op
      char dummy9[sizeof (kim::VariantFilter::StringOperator)];

      // "string value"
      // "info key"
      // string
      char dummy10[sizeof (std::string)];
    };

    /// The size of the largest semantic type.
    enum { size = sizeof (union_type) };

    /// A buffer to store semantic values.
    union
    {
      /// Strongest alignment constraints.
      long double yyalign_me_;
      /// A buffer large enough to store any of the semantic values.
      char yyraw_[size];
    };

    /// Whether the content is built: if defined, the name of the stored type.
    const std::type_info *yytypeid_;
  };

#endif
    /// Backward compatibility (Bison 3.8).
    typedef value_type semantic_type;

    /// Symbol locations.
    typedef location location_type;

    /// Syntax errors thrown from user actions.
    struct syntax_error : std::runtime_error
    {
      syntax_error (const location_type& l, const std::string& m)
        : std::runtime_error (m)
        , location (l)
      {}

      syntax_error (const syntax_error& s)
        : std::runtime_error (s.what ())
        , location (s.location)
      {}

      ~syntax_error () YY_NOEXCEPT YY_NOTHROW;

      location_type location;
    };

    /// Token kinds.
    struct token
    {
      enum token_kind_type
      {
        TOKEN_YYEMPTY = -2,
    TOKEN_YYEOF = 0,               // "end of file"
    TOKEN_YYerror = 1,             // error
    TOKEN_YYUNDEF = 2,             // "invalid token"
    TOKEN_AND = 3,                 // "conjunction operator ('AND' or '&&')"
    TOKEN_OR = 4,                  // "disjunction operator ('or' or '||')"
    TOKEN_NOT = 5,                 // "negation operator ('not' or '!')"
    TOKEN_BLOCK_OPEN = 6,          // "opening parenthesis '('"
    TOKEN_BLOCK_CLOSE = 7,         // "closing parenthesis ')'"
    TOKEN_EQUAL = 8,               // "equality operator ('=')"
    TOKEN_REGEX = 9,               // "regular expression match operator ('~')"
    TOKEN_NOT_EQUAL = 10,          // "inequality operator ('!=' or '<>')"
    TOKEN_LESS_THAN = 11,          // "order relationship operator '<'"
    TOKEN_LESS_OR_EQUAL = 12,      // "order relationship operator '<='"
    TOKEN_GREATER_THAN = 13,       // "order relationship operator '>'"
    TOKEN_GREATER_OR_EQUAL = 14,   // "order relationship operator '>='"
    TOKEN_ON_CHROM = 15,           // "field attribute 'CHROM'"
    TOKEN_ON_POS = 16,             // "field attribute 'POS'"
    TOKEN_ON_ID = 17,              // "field attribute 'ID'"
    TOKEN_ON_QUAL = 18,            // "field attribute 'QUAL'"
    TOKEN_ON_FILTER = 19,          // "field attribute 'FILTER'"
    TOKEN_ON_INFO = 20,            // "field attribute 'INFO'"
    TOKEN_ON_PLOIDY = 21,          // "status 'Ploidy'"
    TOKEN_ON_SNP = 22,             // "status 'SNP' (Single Nucleotide Polymorphism)"
    TOKEN_ON_MULTI_ALLELIC_SNP = 23, // "status 'MSNP' (multi-allelic SNP)"
    TOKEN_ON_SV = 24,              // "status 'SV' (Structural Variant)"
    TOKEN_ON_INDEL = 25,           // "status 'InDel' (Insertion-Deletion)"
    TOKEN_STR_VALUE = 26,          // "string value"
    TOKEN_INFO_FIELD = 27,         // "info key"
    TOKEN_INT_VALUE = 28,          // "integer value"
    TOKEN_NUM_VALUE = 29,          // "number value"
    TOKEN_BOOL_VALUE = 30          // "boolean value"
      };
      /// Backward compatibility alias (Bison 3.6).
      typedef token_kind_type yytokentype;
    };

    /// Token kind, as returned by yylex.
    typedef token::token_kind_type token_kind_type;

    /// Backward compatibility alias (Bison 3.6).
    typedef token_kind_type token_type;

    /// Symbol kinds.
    struct symbol_kind
    {
      enum symbol_kind_type
      {
        YYNTOKENS = 31, ///< Number of tokens.
        S_YYEMPTY = -2,
        S_YYEOF = 0,                             // "end of file"
        S_YYerror = 1,                           // error
        S_YYUNDEF = 2,                           // "invalid token"
        S_AND = 3,                               // "conjunction operator ('AND' or '&&')"
        S_OR = 4,                                // "disjunction operator ('or' or '||')"
        S_NOT = 5,                               // "negation operator ('not' or '!')"
        S_BLOCK_OPEN = 6,                        // "opening parenthesis '('"
        S_BLOCK_CLOSE = 7,                       // "closing parenthesis ')'"
        S_EQUAL = 8,                             // "equality operator ('=')"
        S_REGEX = 9,                             // "regular expression match operator ('~')"
        S_NOT_EQUAL = 10,                        // "inequality operator ('!=' or '<>')"
        S_LESS_THAN = 11,                        // "order relationship operator '<'"
        S_LESS_OR_EQUAL = 12,                    // "order relationship operator '<='"
        S_GREATER_THAN = 13,                     // "order relationship operator '>'"
        S_GREATER_OR_EQUAL = 14,                 // "order relationship operator '>='"
        S_ON_CHROM = 15,                         // "field attribute 'CHROM'"
        S_ON_POS = 16,                           // "field attribute 'POS'"
        S_ON_ID = 17,                            // "field attribute 'ID'"
        S_ON_QUAL = 18,                          // "field attribute 'QUAL'"
        S_ON_FILTER = 19,                        // "field attribute 'FILTER'"
        S_ON_INFO = 20,                          // "field attribute 'INFO'"
        S_ON_PLOIDY = 21,                        // "status 'Ploidy'"
        S_ON_SNP = 22,                           // "status 'SNP' (Single Nucleotide Polymorphism)"
        S_ON_MULTI_ALLELIC_SNP = 23,             // "status 'MSNP' (multi-allelic SNP)"
        S_ON_SV = 24,                            // "status 'SV' (Structural Variant)"
        S_ON_INDEL = 25,                         // "status 'InDel' (Insertion-Deletion)"
        S_STR_VALUE = 26,                        // "string value"
        S_INFO_FIELD = 27,                       // "info key"
        S_INT_VALUE = 28,                        // "integer value"
        S_NUM_VALUE = 29,                        // "number value"
        S_BOOL_VALUE = 30,                       // "boolean value"
        S_YYACCEPT = 31,                         // $accept
        S_expr = 32,                             // expr
        S_filter = 33,                           // filter
        S_simple_expr = 34,                      // simple_expr
        S_any_op = 35,                           // any_op
        S_str_op = 36,                           // str_op
        S_num_op = 37,                           // num_op
        S_bool_op = 38,                          // bool_op
        S_string = 39,                           // string
        S_integer = 40,                          // integer
        S_number = 41,                           // number
        S_boolean = 42                           // boolean
      };
    };

    /// (Internal) symbol kind.
    typedef symbol_kind::symbol_kind_type symbol_kind_type;

    /// The number of tokens.
    static const symbol_kind_type YYNTOKENS = symbol_kind::YYNTOKENS;

    /// A complete symbol.
    ///
    /// Expects its Base type to provide access to the symbol kind
    /// via kind ().
    ///
    /// Provide access to semantic value and location.
    template <typename Base>
    struct basic_symbol : Base
    {
      /// Alias to Base.
      typedef Base super_type;

      /// Default constructor.
      basic_symbol () YY_NOEXCEPT
        : value ()
        , location ()
      {}

#if 201103L <= YY_CPLUSPLUS
      /// Move constructor.
      basic_symbol (basic_symbol&& that)
        : Base (std::move (that))
        , value ()
        , location (std::move (that.location))
      {
        switch (this->kind ())
    {
      case symbol_kind::S_BOOL_VALUE: // "boolean value"
      case symbol_kind::S_filter: // filter
      case symbol_kind::S_simple_expr: // simple_expr
      case symbol_kind::S_boolean: // boolean
        value.move< bool > (std::move (that.value));
        break;

      case symbol_kind::S_NUM_VALUE: // "number value"
      case symbol_kind::S_number: // number
        value.move< double > (std::move (that.value));
        break;

      case symbol_kind::S_INT_VALUE: // "integer value"
      case symbol_kind::S_integer: // integer
        value.move< int > (std::move (that.value));
        break;

      case symbol_kind::S_bool_op: // bool_op
        value.move< kim::VariantFilter::BooleanOperator > (std::move (that.value));
        break;

      case symbol_kind::S_EQUAL: // "equality operator ('=')"
      case symbol_kind::S_REGEX: // "regular expression match operator ('~')"
      case symbol_kind::S_NOT_EQUAL: // "inequality operator ('!=' or '<>')"
      case symbol_kind::S_LESS_THAN: // "order relationship operator '<'"
      case symbol_kind::S_LESS_OR_EQUAL: // "order relationship operator '<='"
      case symbol_kind::S_GREATER_THAN: // "order relationship operator '>'"
      case symbol_kind::S_GREATER_OR_EQUAL: // "order relationship operator '>='"
      case symbol_kind::S_any_op: // any_op
        value.move< kim::VariantFilter::ComparisonOperator > (std::move (that.value));
        break;

      case symbol_kind::S_AND: // "conjunction operator ('AND' or '&&')"
      case symbol_kind::S_OR: // "disjunction operator ('or' or '||')"
      case symbol_kind::S_NOT: // "negation operator ('not' or '!')"
      case symbol_kind::S_BLOCK_OPEN: // "opening parenthesis '('"
      case symbol_kind::S_BLOCK_CLOSE: // "closing parenthesis ')'"
        value.move< kim::VariantFilter::ExpressionOperator > (std::move (that.value));
        break;

      case symbol_kind::S_ON_CHROM: // "field attribute 'CHROM'"
      case symbol_kind::S_ON_POS: // "field attribute 'POS'"
      case symbol_kind::S_ON_ID: // "field attribute 'ID'"
      case symbol_kind::S_ON_QUAL: // "field attribute 'QUAL'"
      case symbol_kind::S_ON_FILTER: // "field attribute 'FILTER'"
      case symbol_kind::S_ON_INFO: // "field attribute 'INFO'"
      case symbol_kind::S_ON_PLOIDY: // "status 'Ploidy'"
      case symbol_kind::S_ON_SNP: // "status 'SNP' (Single Nucleotide Polymorphism)"
      case symbol_kind::S_ON_MULTI_ALLELIC_SNP: // "status 'MSNP' (multi-allelic SNP)"
      case symbol_kind::S_ON_SV: // "status 'SV' (Structural Variant)"
      case symbol_kind::S_ON_INDEL: // "status 'InDel' (Insertion-Deletion)"
        value.move< kim::VariantFilter::FilterOperation > (std::move (that.value));
        break;

      case symbol_kind::S_num_op: // num_op
        value.move< kim::VariantFilter::NumericalOperator > (std::move (that.value));
        break;

      case symbol_kind::S_str_op: // str_op
        value.move< kim::VariantFilter::StringOperator > (std::move (that.value));
        break;

      case symbol_kind::S_STR_VALUE: // "string value"
      case symbol_kind::S_INFO_FIELD: // "info key"
      case symbol_kind::S_string: // string
        value.move< std::string > (std::move (that.value));
        break;

      default:
        break;
    }

      }
#endif

      /// Copy constructor.
      basic_symbol (const basic_symbol& that);

      /// Constructors for typed symbols.
#if 201103L <= YY_CPLUSPLUS
      basic_symbol (typename Base::kind_type t, location_type&& l)
        : Base (t)
        , location (std::move (l))
      {}
#else
      basic_symbol (typename Base::kind_type t, const location_type& l)
        : Base (t)
        , location (l)
      {}
#endif

#if 201103L <= YY_CPLUSPLUS
      basic_symbol (typename Base::kind_type t, bool&& v, location_type&& l)
        : Base (t)
        , value (std::move (v))
        , location (std::move (l))
      {}
#else
      basic_symbol (typename Base::kind_type t, const bool& v, const location_type& l)
        : Base (t)
        , value (v)
        , location (l)
      {}
#endif

#if 201103L <= YY_CPLUSPLUS
      basic_symbol (typename Base::kind_type t, double&& v, location_type&& l)
        : Base (t)
        , value (std::move (v))
        , location (std::move (l))
      {}
#else
      basic_symbol (typename Base::kind_type t, const double& v, const location_type& l)
        : Base (t)
        , value (v)
        , location (l)
      {}
#endif

#if 201103L <= YY_CPLUSPLUS
      basic_symbol (typename Base::kind_type t, int&& v, location_type&& l)
        : Base (t)
        , value (std::move (v))
        , location (std::move (l))
      {}
#else
      basic_symbol (typename Base::kind_type t, const int& v, const location_type& l)
        : Base (t)
        , value (v)
        , location (l)
      {}
#endif

#if 201103L <= YY_CPLUSPLUS
      basic_symbol (typename Base::kind_type t, kim::VariantFilter::BooleanOperator&& v, location_type&& l)
        : Base (t)
        , value (std::move (v))
        , location (std::move (l))
      {}
#else
      basic_symbol (typename Base::kind_type t, const kim::VariantFilter::BooleanOperator& v, const location_type& l)
        : Base (t)
        , value (v)
        , location (l)
      {}
#endif

#if 201103L <= YY_CPLUSPLUS
      basic_symbol (typename Base::kind_type t, kim::VariantFilter::ComparisonOperator&& v, location_type&& l)
        : Base (t)
        , value (std::move (v))
        , location (std::move (l))
      {}
#else
      basic_symbol (typename Base::kind_type t, const kim::VariantFilter::ComparisonOperator& v, const location_type& l)
        : Base (t)
        , value (v)
        , location (l)
      {}
#endif

#if 201103L <= YY_CPLUSPLUS
      basic_symbol (typename Base::kind_type t, kim::VariantFilter::ExpressionOperator&& v, location_type&& l)
        : Base (t)
        , value (std::move (v))
        , location (std::move (l))
      {}
#else
      basic_symbol (typename Base::kind_type t, const kim::VariantFilter::ExpressionOperator& v, const location_type& l)
        : Base (t)
        , value (v)
        , location (l)
      {}
#endif

#if 201103L <= YY_CPLUSPLUS
      basic_symbol (typename Base::kind_type t, kim::VariantFilter::FilterOperation&& v, location_type&& l)
        : Base (t)
        , value (std::move (v))
        , location (std::move (l))
      {}
#else
      basic_symbol (typename Base::kind_type t, const kim::VariantFilter::FilterOperation& v, const location_type& l)
        : Base (t)
        , value (v)
        , location (l)
      {}
#endif

#if 201103L <= YY_CPLUSPLUS
      basic_symbol (typename Base::kind_type t, kim::VariantFilter::NumericalOperator&& v, location_type&& l)
        : Base (t)
        , value (std::move (v))
        , location (std::move (l))
      {}
#else
      basic_symbol (typename Base::kind_type t, const kim::VariantFilter::NumericalOperator& v, const location_type& l)
        : Base (t)
        , value (v)
        , location (l)
      {}
#endif

#if 201103L <= YY_CPLUSPLUS
      basic_symbol (typename Base::kind_type t, kim::VariantFilter::StringOperator&& v, location_type&& l)
        : Base (t)
        , value (std::move (v))
        , location (std::move (l))
      {}
#else
      basic_symbol (typename Base::kind_type t, const kim::VariantFilter::StringOperator& v, const location_type& l)
        : Base (t)
        , value (v)
        , location (l)
      {}
#endif

#if 201103L <= YY_CPLUSPLUS
      basic_symbol (typename Base::kind_type t, std::string&& v, location_type&& l)
        : Base (t)
        , value (std::move (v))
        , location (std::move (l))
      {}
#else
      basic_symbol (typename Base::kind_type t, const std::string& v, const location_type& l)
        : Base (t)
        , value (v)
        , location (l)
      {}
#endif

      /// Destroy the symbol.
      ~basic_symbol ()
      {
        clear ();
      }



      /// Destroy contents, and record that is empty.
      void clear () YY_NOEXCEPT
      {
        // User destructor.
        symbol_kind_type yykind = this->kind ();
        basic_symbol<Base>& yysym = *this;
        (void) yysym;
        switch (yykind)
        {
       default:
          break;
        }

        // Value type destructor.
switch (yykind)
    {
      case symbol_kind::S_BOOL_VALUE: // "boolean value"
      case symbol_kind::S_filter: // filter
      case symbol_kind::S_simple_expr: // simple_expr
      case symbol_kind::S_boolean: // boolean
        value.template destroy< bool > ();
        break;

      case symbol_kind::S_NUM_VALUE: // "number value"
      case symbol_kind::S_number: // number
        value.template destroy< double > ();
        break;

      case symbol_kind::S_INT_VALUE: // "integer value"
      case symbol_kind::S_integer: // integer
        value.template destroy< int > ();
        break;

      case symbol_kind::S_bool_op: // bool_op
        value.template destroy< kim::VariantFilter::BooleanOperator > ();
        break;

      case symbol_kind::S_EQUAL: // "equality operator ('=')"
      case symbol_kind::S_REGEX: // "regular expression match operator ('~')"
      case symbol_kind::S_NOT_EQUAL: // "inequality operator ('!=' or '<>')"
      case symbol_kind::S_LESS_THAN: // "order relationship operator '<'"
      case symbol_kind::S_LESS_OR_EQUAL: // "order relationship operator '<='"
      case symbol_kind::S_GREATER_THAN: // "order relationship operator '>'"
      case symbol_kind::S_GREATER_OR_EQUAL: // "order relationship operator '>='"
      case symbol_kind::S_any_op: // any_op
        value.template destroy< kim::VariantFilter::ComparisonOperator > ();
        break;

      case symbol_kind::S_AND: // "conjunction operator ('AND' or '&&')"
      case symbol_kind::S_OR: // "disjunction operator ('or' or '||')"
      case symbol_kind::S_NOT: // "negation operator ('not' or '!')"
      case symbol_kind::S_BLOCK_OPEN: // "opening parenthesis '('"
      case symbol_kind::S_BLOCK_CLOSE: // "closing parenthesis ')'"
        value.template destroy< kim::VariantFilter::ExpressionOperator > ();
        break;

      case symbol_kind::S_ON_CHROM: // "field attribute 'CHROM'"
      case symbol_kind::S_ON_POS: // "field attribute 'POS'"
      case symbol_kind::S_ON_ID: // "field attribute 'ID'"
      case symbol_kind::S_ON_QUAL: // "field attribute 'QUAL'"
      case symbol_kind::S_ON_FILTER: // "field attribute 'FILTER'"
      case symbol_kind::S_ON_INFO: // "field attribute 'INFO'"
      case symbol_kind::S_ON_PLOIDY: // "status 'Ploidy'"
      case symbol_kind::S_ON_SNP: // "status 'SNP' (Single Nucleotide Polymorphism)"
      case symbol_kind::S_ON_MULTI_ALLELIC_SNP: // "status 'MSNP' (multi-allelic SNP)"
      case symbol_kind::S_ON_SV: // "status 'SV' (Structural Variant)"
      case symbol_kind::S_ON_INDEL: // "status 'InDel' (Insertion-Deletion)"
        value.template destroy< kim::VariantFilter::FilterOperation > ();
        break;

      case symbol_kind::S_num_op: // num_op
        value.template destroy< kim::VariantFilter::NumericalOperator > ();
        break;

      case symbol_kind::S_str_op: // str_op
        value.template destroy< kim::VariantFilter::StringOperator > ();
        break;

      case symbol_kind::S_STR_VALUE: // "string value"
      case symbol_kind::S_INFO_FIELD: // "info key"
      case symbol_kind::S_string: // string
        value.template destroy< std::string > ();
        break;

      default:
        break;
    }

        Base::clear ();
      }

      /// The user-facing name of this symbol.
      const char *name () const YY_NOEXCEPT
      {
        return VariantFilterParser::symbol_name (this->kind ());
      }

      /// Backward compatibility (Bison 3.6).
      symbol_kind_type type_get () const YY_NOEXCEPT;

      /// Whether empty.
      bool empty () const YY_NOEXCEPT;

      /// Destructive move, \a s is emptied into this.
      void move (basic_symbol& s);

      /// The semantic value.
      value_type value;

      /// The location.
      location_type location;

    private:
#if YY_CPLUSPLUS < 201103L
      /// Assignment operator.
      basic_symbol& operator= (const basic_symbol& that);
#endif
    };

    /// Type access provider for token (enum) based symbols.
    struct by_kind
    {
      /// The symbol kind as needed by the constructor.
      typedef token_kind_type kind_type;

      /// Default constructor.
      by_kind () YY_NOEXCEPT;

#if 201103L <= YY_CPLUSPLUS
      /// Move constructor.
      by_kind (by_kind&& that) YY_NOEXCEPT;
#endif

      /// Copy constructor.
      by_kind (const by_kind& that) YY_NOEXCEPT;

      /// Constructor from (external) token numbers.
      by_kind (kind_type t) YY_NOEXCEPT;



      /// Record that this symbol is empty.
      void clear () YY_NOEXCEPT;

      /// Steal the symbol kind from \a that.
      void move (by_kind& that);

      /// The (internal) type number (corresponding to \a type).
      /// \a empty when empty.
      symbol_kind_type kind () const YY_NOEXCEPT;

      /// Backward compatibility (Bison 3.6).
      symbol_kind_type type_get () const YY_NOEXCEPT;

      /// The symbol kind.
      /// \a S_YYEMPTY when empty.
      symbol_kind_type kind_;
    };

    /// Backward compatibility for a private implementation detail (Bison 3.6).
    typedef by_kind by_type;

    /// "External" symbols: returned by the scanner.
    struct symbol_type : basic_symbol<by_kind>
    {
      /// Superclass.
      typedef basic_symbol<by_kind> super_type;

      /// Empty symbol.
      symbol_type () YY_NOEXCEPT {}

      /// Constructor for valueless symbols, and symbols from each type.
#if 201103L <= YY_CPLUSPLUS
      symbol_type (int tok, location_type l)
        : super_type (token_kind_type (tok), std::move (l))
#else
      symbol_type (int tok, const location_type& l)
        : super_type (token_kind_type (tok), l)
#endif
      {
#if !defined _MSC_VER || defined __clang__
        YY_ASSERT (tok == token::TOKEN_YYEOF
                   || (token::TOKEN_YYerror <= tok && tok <= token::TOKEN_YYUNDEF));
#endif
      }
#if 201103L <= YY_CPLUSPLUS
      symbol_type (int tok, bool v, location_type l)
        : super_type (token_kind_type (tok), std::move (v), std::move (l))
#else
      symbol_type (int tok, const bool& v, const location_type& l)
        : super_type (token_kind_type (tok), v, l)
#endif
      {
#if !defined _MSC_VER || defined __clang__
        YY_ASSERT (tok == token::TOKEN_BOOL_VALUE);
#endif
      }
#if 201103L <= YY_CPLUSPLUS
      symbol_type (int tok, double v, location_type l)
        : super_type (token_kind_type (tok), std::move (v), std::move (l))
#else
      symbol_type (int tok, const double& v, const location_type& l)
        : super_type (token_kind_type (tok), v, l)
#endif
      {
#if !defined _MSC_VER || defined __clang__
        YY_ASSERT (tok == token::TOKEN_NUM_VALUE);
#endif
      }
#if 201103L <= YY_CPLUSPLUS
      symbol_type (int tok, int v, location_type l)
        : super_type (token_kind_type (tok), std::move (v), std::move (l))
#else
      symbol_type (int tok, const int& v, const location_type& l)
        : super_type (token_kind_type (tok), v, l)
#endif
      {
#if !defined _MSC_VER || defined __clang__
        YY_ASSERT (tok == token::TOKEN_INT_VALUE);
#endif
      }
#if 201103L <= YY_CPLUSPLUS
      symbol_type (int tok, kim::VariantFilter::ComparisonOperator v, location_type l)
        : super_type (token_kind_type (tok), std::move (v), std::move (l))
#else
      symbol_type (int tok, const kim::VariantFilter::ComparisonOperator& v, const location_type& l)
        : super_type (token_kind_type (tok), v, l)
#endif
      {
#if !defined _MSC_VER || defined __clang__
        YY_ASSERT ((token::TOKEN_EQUAL <= tok && tok <= token::TOKEN_GREATER_OR_EQUAL));
#endif
      }
#if 201103L <= YY_CPLUSPLUS
      symbol_type (int tok, kim::VariantFilter::ExpressionOperator v, location_type l)
        : super_type (token_kind_type (tok), std::move (v), std::move (l))
#else
      symbol_type (int tok, const kim::VariantFilter::ExpressionOperator& v, const location_type& l)
        : super_type (token_kind_type (tok), v, l)
#endif
      {
#if !defined _MSC_VER || defined __clang__
        YY_ASSERT ((token::TOKEN_AND <= tok && tok <= token::TOKEN_BLOCK_CLOSE));
#endif
      }
#if 201103L <= YY_CPLUSPLUS
      symbol_type (int tok, kim::VariantFilter::FilterOperation v, location_type l)
        : super_type (token_kind_type (tok), std::move (v), std::move (l))
#else
      symbol_type (int tok, const kim::VariantFilter::FilterOperation& v, const location_type& l)
        : super_type (token_kind_type (tok), v, l)
#endif
      {
#if !defined _MSC_VER || defined __clang__
        YY_ASSERT ((token::TOKEN_ON_CHROM <= tok && tok <= token::TOKEN_ON_INDEL));
#endif
      }
#if 201103L <= YY_CPLUSPLUS
      symbol_type (int tok, std::string v, location_type l)
        : super_type (token_kind_type (tok), std::move (v), std::move (l))
#else
      symbol_type (int tok, const std::string& v, const location_type& l)
        : super_type (token_kind_type (tok), v, l)
#endif
      {
#if !defined _MSC_VER || defined __clang__
        YY_ASSERT ((token::TOKEN_STR_VALUE <= tok && tok <= token::TOKEN_INFO_FIELD));
#endif
      }
    };

    /// Build a parser object.
    VariantFilterParser (VariantFilterDriver &driver_yyarg);
    virtual ~VariantFilterParser ();

#if 201103L <= YY_CPLUSPLUS
    /// Non copyable.
    VariantFilterParser (const VariantFilterParser&) = delete;
    /// Non copyable.
    VariantFilterParser& operator= (const VariantFilterParser&) = delete;
#endif

    /// Parse.  An alias for parse ().
    /// \returns  0 iff parsing succeeded.
    int operator() ();

    /// Parse.
    /// \returns  0 iff parsing succeeded.
    virtual int parse ();

#if YYDEBUG
    /// The current debugging stream.
    std::ostream& debug_stream () const YY_ATTRIBUTE_PURE;
    /// Set the current debugging stream.
    void set_debug_stream (std::ostream &);

    /// Type for debugging levels.
    typedef int debug_level_type;
    /// The current debugging level.
    debug_level_type debug_level () const YY_ATTRIBUTE_PURE;
    /// Set the current debugging level.
    void set_debug_level (debug_level_type l);
#endif

    /// Report a syntax error.
    /// \param loc    where the syntax error is found.
    /// \param msg    a description of the syntax error.
    virtual void error (const location_type& loc, const std::string& msg);

    /// Report a syntax error.
    void error (const syntax_error& err);

    /// The user-facing name of the symbol whose (internal) number is
    /// YYSYMBOL.  No bounds checking.
    static const char *symbol_name (symbol_kind_type yysymbol);

    // Implementation of make_symbol for each token kind.
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_YYEOF (location_type l)
      {
        return symbol_type (token::TOKEN_YYEOF, std::move (l));
      }
#else
      static
      symbol_type
      make_YYEOF (const location_type& l)
      {
        return symbol_type (token::TOKEN_YYEOF, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_YYerror (location_type l)
      {
        return symbol_type (token::TOKEN_YYerror, std::move (l));
      }
#else
      static
      symbol_type
      make_YYerror (const location_type& l)
      {
        return symbol_type (token::TOKEN_YYerror, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_YYUNDEF (location_type l)
      {
        return symbol_type (token::TOKEN_YYUNDEF, std::move (l));
      }
#else
      static
      symbol_type
      make_YYUNDEF (const location_type& l)
      {
        return symbol_type (token::TOKEN_YYUNDEF, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_AND (kim::VariantFilter::ExpressionOperator v, location_type l)
      {
        return symbol_type (token::TOKEN_AND, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_AND (const kim::VariantFilter::ExpressionOperator& v, const location_type& l)
      {
        return symbol_type (token::TOKEN_AND, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_OR (kim::VariantFilter::ExpressionOperator v, location_type l)
      {
        return symbol_type (token::TOKEN_OR, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_OR (const kim::VariantFilter::ExpressionOperator& v, const location_type& l)
      {
        return symbol_type (token::TOKEN_OR, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_NOT (kim::VariantFilter::ExpressionOperator v, location_type l)
      {
        return symbol_type (token::TOKEN_NOT, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_NOT (const kim::VariantFilter::ExpressionOperator& v, const location_type& l)
      {
        return symbol_type (token::TOKEN_NOT, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_BLOCK_OPEN (kim::VariantFilter::ExpressionOperator v, location_type l)
      {
        return symbol_type (token::TOKEN_BLOCK_OPEN, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_BLOCK_OPEN (const kim::VariantFilter::ExpressionOperator& v, const location_type& l)
      {
        return symbol_type (token::TOKEN_BLOCK_OPEN, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_BLOCK_CLOSE (kim::VariantFilter::ExpressionOperator v, location_type l)
      {
        return symbol_type (token::TOKEN_BLOCK_CLOSE, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_BLOCK_CLOSE (const kim::VariantFilter::ExpressionOperator& v, const location_type& l)
      {
        return symbol_type (token::TOKEN_BLOCK_CLOSE, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_EQUAL (kim::VariantFilter::ComparisonOperator v, location_type l)
      {
        return symbol_type (token::TOKEN_EQUAL, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_EQUAL (const kim::VariantFilter::ComparisonOperator& v, const location_type& l)
      {
        return symbol_type (token::TOKEN_EQUAL, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_REGEX (kim::VariantFilter::ComparisonOperator v, location_type l)
      {
        return symbol_type (token::TOKEN_REGEX, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_REGEX (const kim::VariantFilter::ComparisonOperator& v, const location_type& l)
      {
        return symbol_type (token::TOKEN_REGEX, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_NOT_EQUAL (kim::VariantFilter::ComparisonOperator v, location_type l)
      {
        return symbol_type (token::TOKEN_NOT_EQUAL, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_NOT_EQUAL (const kim::VariantFilter::ComparisonOperator& v, const location_type& l)
      {
        return symbol_type (token::TOKEN_NOT_EQUAL, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_LESS_THAN (kim::VariantFilter::ComparisonOperator v, location_type l)
      {
        return symbol_type (token::TOKEN_LESS_THAN, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_LESS_THAN (const kim::VariantFilter::ComparisonOperator& v, const location_type& l)
      {
        return symbol_type (token::TOKEN_LESS_THAN, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_LESS_OR_EQUAL (kim::VariantFilter::ComparisonOperator v, location_type l)
      {
        return symbol_type (token::TOKEN_LESS_OR_EQUAL, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_LESS_OR_EQUAL (const kim::VariantFilter::ComparisonOperator& v, const location_type& l)
      {
        return symbol_type (token::TOKEN_LESS_OR_EQUAL, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_GREATER_THAN (kim::VariantFilter::ComparisonOperator v, location_type l)
      {
        return symbol_type (token::TOKEN_GREATER_THAN, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_GREATER_THAN (const kim::VariantFilter::ComparisonOperator& v, const location_type& l)
      {
        return symbol_type (token::TOKEN_GREATER_THAN, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_GREATER_OR_EQUAL (kim::VariantFilter::ComparisonOperator v, location_type l)
      {
        return symbol_type (token::TOKEN_GREATER_OR_EQUAL, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_GREATER_OR_EQUAL (const kim::VariantFilter::ComparisonOperator& v, const location_type& l)
      {
        return symbol_type (token::TOKEN_GREATER_OR_EQUAL, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_ON_CHROM (kim::VariantFilter::FilterOperation v, location_type l)
      {
        return symbol_type (token::TOKEN_ON_CHROM, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_ON_CHROM (const kim::VariantFilter::FilterOperation& v, const location_type& l)
      {
        return symbol_type (token::TOKEN_ON_CHROM, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_ON_POS (kim::VariantFilter::FilterOperation v, location_type l)
      {
        return symbol_type (token::TOKEN_ON_POS, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_ON_POS (const kim::VariantFilter::FilterOperation& v, const location_type& l)
      {
        return symbol_type (token::TOKEN_ON_POS, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_ON_ID (kim::VariantFilter::FilterOperation v, location_type l)
      {
        return symbol_type (token::TOKEN_ON_ID, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_ON_ID (const kim::VariantFilter::FilterOperation& v, const location_type& l)
      {
        return symbol_type (token::TOKEN_ON_ID, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_ON_QUAL (kim::VariantFilter::FilterOperation v, location_type l)
      {
        return symbol_type (token::TOKEN_ON_QUAL, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_ON_QUAL (const kim::VariantFilter::FilterOperation& v, const location_type& l)
      {
        return symbol_type (token::TOKEN_ON_QUAL, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_ON_FILTER (kim::VariantFilter::FilterOperation v, location_type l)
      {
        return symbol_type (token::TOKEN_ON_FILTER, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_ON_FILTER (const kim::VariantFilter::FilterOperation& v, const location_type& l)
      {
        return symbol_type (token::TOKEN_ON_FILTER, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_ON_INFO (kim::VariantFilter::FilterOperation v, location_type l)
      {
        return symbol_type (token::TOKEN_ON_INFO, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_ON_INFO (const kim::VariantFilter::FilterOperation& v, const location_type& l)
      {
        return symbol_type (token::TOKEN_ON_INFO, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_ON_PLOIDY (kim::VariantFilter::FilterOperation v, location_type l)
      {
        return symbol_type (token::TOKEN_ON_PLOIDY, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_ON_PLOIDY (const kim::VariantFilter::FilterOperation& v, const location_type& l)
      {
        return symbol_type (token::TOKEN_ON_PLOIDY, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_ON_SNP (kim::VariantFilter::FilterOperation v, location_type l)
      {
        return symbol_type (token::TOKEN_ON_SNP, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_ON_SNP (const kim::VariantFilter::FilterOperation& v, const location_type& l)
      {
        return symbol_type (token::TOKEN_ON_SNP, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_ON_MULTI_ALLELIC_SNP (kim::VariantFilter::FilterOperation v, location_type l)
      {
        return symbol_type (token::TOKEN_ON_MULTI_ALLELIC_SNP, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_ON_MULTI_ALLELIC_SNP (const kim::VariantFilter::FilterOperation& v, const location_type& l)
      {
        return symbol_type (token::TOKEN_ON_MULTI_ALLELIC_SNP, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_ON_SV (kim::VariantFilter::FilterOperation v, location_type l)
      {
        return symbol_type (token::TOKEN_ON_SV, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_ON_SV (const kim::VariantFilter::FilterOperation& v, const location_type& l)
      {
        return symbol_type (token::TOKEN_ON_SV, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_ON_INDEL (kim::VariantFilter::FilterOperation v, location_type l)
      {
        return symbol_type (token::TOKEN_ON_INDEL, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_ON_INDEL (const kim::VariantFilter::FilterOperation& v, const location_type& l)
      {
        return symbol_type (token::TOKEN_ON_INDEL, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_STR_VALUE (std::string v, location_type l)
      {
        return symbol_type (token::TOKEN_STR_VALUE, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_STR_VALUE (const std::string& v, const location_type& l)
      {
        return symbol_type (token::TOKEN_STR_VALUE, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_INFO_FIELD (std::string v, location_type l)
      {
        return symbol_type (token::TOKEN_INFO_FIELD, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_INFO_FIELD (const std::string& v, const location_type& l)
      {
        return symbol_type (token::TOKEN_INFO_FIELD, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_INT_VALUE (int v, location_type l)
      {
        return symbol_type (token::TOKEN_INT_VALUE, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_INT_VALUE (const int& v, const location_type& l)
      {
        return symbol_type (token::TOKEN_INT_VALUE, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_NUM_VALUE (double v, location_type l)
      {
        return symbol_type (token::TOKEN_NUM_VALUE, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_NUM_VALUE (const double& v, const location_type& l)
      {
        return symbol_type (token::TOKEN_NUM_VALUE, v, l);
      }
#endif
#if 201103L <= YY_CPLUSPLUS
      static
      symbol_type
      make_BOOL_VALUE (bool v, location_type l)
      {
        return symbol_type (token::TOKEN_BOOL_VALUE, std::move (v), std::move (l));
      }
#else
      static
      symbol_type
      make_BOOL_VALUE (const bool& v, const location_type& l)
      {
        return symbol_type (token::TOKEN_BOOL_VALUE, v, l);
      }
#endif


    class context
    {
    public:
      context (const VariantFilterParser& yyparser, const symbol_type& yyla);
      const symbol_type& lookahead () const YY_NOEXCEPT { return yyla_; }
      symbol_kind_type token () const YY_NOEXCEPT { return yyla_.kind (); }
      const location_type& location () const YY_NOEXCEPT { return yyla_.location; }

      /// Put in YYARG at most YYARGN of the expected tokens, and return the
      /// number of tokens stored in YYARG.  If YYARG is null, return the
      /// number of expected tokens (guaranteed to be less than YYNTOKENS).
      int expected_tokens (symbol_kind_type yyarg[], int yyargn) const;

    private:
      const VariantFilterParser& yyparser_;
      const symbol_type& yyla_;
    };

  private:
#if YY_CPLUSPLUS < 201103L
    /// Non copyable.
    VariantFilterParser (const VariantFilterParser&);
    /// Non copyable.
    VariantFilterParser& operator= (const VariantFilterParser&);
#endif

    /// Check the lookahead yytoken.
    /// \returns  true iff the token will be eventually shifted.
    bool yy_lac_check_ (symbol_kind_type yytoken) const;
    /// Establish the initial context if no initial context currently exists.
    /// \returns  true iff the token will be eventually shifted.
    bool yy_lac_establish_ (symbol_kind_type yytoken);
    /// Discard any previous initial lookahead context because of event.
    /// \param event  the event which caused the lookahead to be discarded.
    ///               Only used for debbuging output.
    void yy_lac_discard_ (const char* event);

    /// Stored state numbers (used for stacks).
    typedef signed char state_type;

    /// The arguments of the error message.
    int yy_syntax_error_arguments_ (const context& yyctx,
                                    symbol_kind_type yyarg[], int yyargn) const;

    /// Generate an error message.
    /// \param yyctx     the context in which the error occurred.
    virtual std::string yysyntax_error_ (const context& yyctx) const;
    /// Compute post-reduction state.
    /// \param yystate   the current state
    /// \param yysym     the nonterminal to push on the stack
    static state_type yy_lr_goto_state_ (state_type yystate, int yysym);

    /// Whether the given \c yypact_ value indicates a defaulted state.
    /// \param yyvalue   the value to check
    static bool yy_pact_value_is_default_ (int yyvalue) YY_NOEXCEPT;

    /// Whether the given \c yytable_ value indicates a syntax error.
    /// \param yyvalue   the value to check
    static bool yy_table_value_is_error_ (int yyvalue) YY_NOEXCEPT;

    static const signed char yypact_ninf_;
    static const signed char yytable_ninf_;

    /// Convert a scanner token kind \a t to a symbol kind.
    /// In theory \a t should be a token_kind_type, but character literals
    /// are valid, yet not members of the token_kind_type enum.
    static symbol_kind_type yytranslate_ (int t) YY_NOEXCEPT;



    // Tables.
    // YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
    // STATE-NUM.
    static const signed char yypact_[];

    // YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
    // Performed when YYTABLE does not specify something else to do.  Zero
    // means the default is an error.
    static const signed char yydefact_[];

    // YYPGOTO[NTERM-NUM].
    static const signed char yypgoto_[];

    // YYDEFGOTO[NTERM-NUM].
    static const signed char yydefgoto_[];

    // YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
    // positive, shift that token.  If negative, reduce the rule whose
    // number is the opposite.  If YYTABLE_NINF, syntax error.
    static const signed char yytable_[];

    static const signed char yycheck_[];

    // YYSTOS[STATE-NUM] -- The symbol kind of the accessing symbol of
    // state STATE-NUM.
    static const signed char yystos_[];

    // YYR1[RULE-NUM] -- Symbol kind of the left-hand side of rule RULE-NUM.
    static const signed char yyr1_[];

    // YYR2[RULE-NUM] -- Number of symbols on the right-hand side of rule RULE-NUM.
    static const signed char yyr2_[];


#if YYDEBUG
    // YYRLINE[YYN] -- Source line where rule number YYN was defined.
    static const unsigned char yyrline_[];
    /// Report on the debug stream that the rule \a r is going to be reduced.
    virtual void yy_reduce_print_ (int r) const;
    /// Print the state stack on the debug stream.
    virtual void yy_stack_print_ () const;

    /// Debugging level.
    int yydebug_;
    /// Debug stream.
    std::ostream* yycdebug_;

    /// \brief Display a symbol kind, value and location.
    /// \param yyo    The output stream.
    /// \param yysym  The symbol.
    template <typename Base>
    void yy_print_ (std::ostream& yyo, const basic_symbol<Base>& yysym) const;
#endif

    /// \brief Reclaim the memory associated to a symbol.
    /// \param yymsg     Why this token is reclaimed.
    ///                  If null, print nothing.
    /// \param yysym     The symbol.
    template <typename Base>
    void yy_destroy_ (const char* yymsg, basic_symbol<Base>& yysym) const;

  private:
    /// Type access provider for state based symbols.
    struct by_state
    {
      /// Default constructor.
      by_state () YY_NOEXCEPT;

      /// The symbol kind as needed by the constructor.
      typedef state_type kind_type;

      /// Constructor.
      by_state (kind_type s) YY_NOEXCEPT;

      /// Copy constructor.
      by_state (const by_state& that) YY_NOEXCEPT;

      /// Record that this symbol is empty.
      void clear () YY_NOEXCEPT;

      /// Steal the symbol kind from \a that.
      void move (by_state& that);

      /// The symbol kind (corresponding to \a state).
      /// \a symbol_kind::S_YYEMPTY when empty.
      symbol_kind_type kind () const YY_NOEXCEPT;

      /// The state number used to denote an empty symbol.
      /// We use the initial state, as it does not have a value.
      enum { empty_state = 0 };

      /// The state.
      /// \a empty when empty.
      state_type state;
    };

    /// "Internal" symbol: element of the stack.
    struct stack_symbol_type : basic_symbol<by_state>
    {
      /// Superclass.
      typedef basic_symbol<by_state> super_type;
      /// Construct an empty symbol.
      stack_symbol_type ();
      /// Move or copy construction.
      stack_symbol_type (YY_RVREF (stack_symbol_type) that);
      /// Steal the contents from \a sym to build this.
      stack_symbol_type (state_type s, YY_MOVE_REF (symbol_type) sym);
#if YY_CPLUSPLUS < 201103L
      /// Assignment, needed by push_back by some old implementations.
      /// Moves the contents of that.
      stack_symbol_type& operator= (stack_symbol_type& that);

      /// Assignment, needed by push_back by other implementations.
      /// Needed by some other old implementations.
      stack_symbol_type& operator= (const stack_symbol_type& that);
#endif
    };

    /// A stack with random access from its top.
    template <typename T, typename S = std::vector<T> >
    class stack
    {
    public:
      // Hide our reversed order.
      typedef typename S::iterator iterator;
      typedef typename S::const_iterator const_iterator;
      typedef typename S::size_type size_type;
      typedef typename std::ptrdiff_t index_type;

      stack (size_type n = 200) YY_NOEXCEPT
        : seq_ (n)
      {}

#if 201103L <= YY_CPLUSPLUS
      /// Non copyable.
      stack (const stack&) = delete;
      /// Non copyable.
      stack& operator= (const stack&) = delete;
#endif

      /// Random access.
      ///
      /// Index 0 returns the topmost element.
      const T&
      operator[] (index_type i) const
      {
        return seq_[size_type (size () - 1 - i)];
      }

      /// Random access.
      ///
      /// Index 0 returns the topmost element.
      T&
      operator[] (index_type i)
      {
        return seq_[size_type (size () - 1 - i)];
      }

      /// Steal the contents of \a t.
      ///
      /// Close to move-semantics.
      void
      push (YY_MOVE_REF (T) t)
      {
        seq_.push_back (T ());
        operator[] (0).move (t);
      }

      /// Pop elements from the stack.
      void
      pop (std::ptrdiff_t n = 1) YY_NOEXCEPT
      {
        for (; 0 < n; --n)
          seq_.pop_back ();
      }

      /// Pop all elements from the stack.
      void
      clear () YY_NOEXCEPT
      {
        seq_.clear ();
      }

      /// Number of elements on the stack.
      index_type
      size () const YY_NOEXCEPT
      {
        return index_type (seq_.size ());
      }

      /// Iterator on top of the stack (going downwards).
      const_iterator
      begin () const YY_NOEXCEPT
      {
        return seq_.begin ();
      }

      /// Bottom of the stack.
      const_iterator
      end () const YY_NOEXCEPT
      {
        return seq_.end ();
      }

      /// Present a slice of the top of a stack.
      class slice
      {
      public:
        slice (const stack& stack, index_type range) YY_NOEXCEPT
          : stack_ (stack)
          , range_ (range)
        {}

        const T&
        operator[] (index_type i) const
        {
          return stack_[range_ - i];
        }

      private:
        const stack& stack_;
        index_type range_;
      };

    private:
#if YY_CPLUSPLUS < 201103L
      /// Non copyable.
      stack (const stack&);
      /// Non copyable.
      stack& operator= (const stack&);
#endif
      /// The wrapped container.
      S seq_;
    };


    /// Stack type.
    typedef stack<stack_symbol_type> stack_type;

    /// The stack.
    stack_type yystack_;
    /// The stack for LAC.
    /// Logically, the yy_lac_stack's lifetime is confined to the function
    /// yy_lac_check_. We just store it as a member of this class to hold
    /// on to the memory and to avoid frequent reallocations.
    /// Since yy_lac_check_ is const, this member must be mutable.
    mutable std::vector<state_type> yylac_stack_;
    /// Whether an initial LAC context was established.
    bool yy_lac_established_;


    /// Push a new state on the stack.
    /// \param m    a debug message to display
    ///             if null, no trace is output.
    /// \param sym  the symbol
    /// \warning the contents of \a s.value is stolen.
    void yypush_ (const char* m, YY_MOVE_REF (stack_symbol_type) sym);

    /// Push a new look ahead token on the state on the stack.
    /// \param m    a debug message to display
    ///             if null, no trace is output.
    /// \param s    the state
    /// \param sym  the symbol (for its value and location).
    /// \warning the contents of \a sym.value is stolen.
    void yypush_ (const char* m, state_type s, YY_MOVE_REF (symbol_type) sym);

    /// Pop \a n symbols from the stack.
    void yypop_ (int n = 1) YY_NOEXCEPT;

    /// Constants.
    enum
    {
      yylast_ = 67,     ///< Last index in yytable_.
      yynnts_ = 12,  ///< Number of nonterminal symbols.
      yyfinal_ = 41 ///< Termination state number.
    };


    // User arguments.
    VariantFilterDriver &driver;

  };

  inline
  VariantFilterParser::symbol_kind_type
  VariantFilterParser::yytranslate_ (int t) YY_NOEXCEPT
  {
    return static_cast<symbol_kind_type> (t);
  }

  // basic_symbol.
  template <typename Base>
  VariantFilterParser::basic_symbol<Base>::basic_symbol (const basic_symbol& that)
    : Base (that)
    , value ()
    , location (that.location)
  {
    switch (this->kind ())
    {
      case symbol_kind::S_BOOL_VALUE: // "boolean value"
      case symbol_kind::S_filter: // filter
      case symbol_kind::S_simple_expr: // simple_expr
      case symbol_kind::S_boolean: // boolean
        value.copy< bool > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_NUM_VALUE: // "number value"
      case symbol_kind::S_number: // number
        value.copy< double > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_INT_VALUE: // "integer value"
      case symbol_kind::S_integer: // integer
        value.copy< int > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_bool_op: // bool_op
        value.copy< kim::VariantFilter::BooleanOperator > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_EQUAL: // "equality operator ('=')"
      case symbol_kind::S_REGEX: // "regular expression match operator ('~')"
      case symbol_kind::S_NOT_EQUAL: // "inequality operator ('!=' or '<>')"
      case symbol_kind::S_LESS_THAN: // "order relationship operator '<'"
      case symbol_kind::S_LESS_OR_EQUAL: // "order relationship operator '<='"
      case symbol_kind::S_GREATER_THAN: // "order relationship operator '>'"
      case symbol_kind::S_GREATER_OR_EQUAL: // "order relationship operator '>='"
      case symbol_kind::S_any_op: // any_op
        value.copy< kim::VariantFilter::ComparisonOperator > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_AND: // "conjunction operator ('AND' or '&&')"
      case symbol_kind::S_OR: // "disjunction operator ('or' or '||')"
      case symbol_kind::S_NOT: // "negation operator ('not' or '!')"
      case symbol_kind::S_BLOCK_OPEN: // "opening parenthesis '('"
      case symbol_kind::S_BLOCK_CLOSE: // "closing parenthesis ')'"
        value.copy< kim::VariantFilter::ExpressionOperator > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_ON_CHROM: // "field attribute 'CHROM'"
      case symbol_kind::S_ON_POS: // "field attribute 'POS'"
      case symbol_kind::S_ON_ID: // "field attribute 'ID'"
      case symbol_kind::S_ON_QUAL: // "field attribute 'QUAL'"
      case symbol_kind::S_ON_FILTER: // "field attribute 'FILTER'"
      case symbol_kind::S_ON_INFO: // "field attribute 'INFO'"
      case symbol_kind::S_ON_PLOIDY: // "status 'Ploidy'"
      case symbol_kind::S_ON_SNP: // "status 'SNP' (Single Nucleotide Polymorphism)"
      case symbol_kind::S_ON_MULTI_ALLELIC_SNP: // "status 'MSNP' (multi-allelic SNP)"
      case symbol_kind::S_ON_SV: // "status 'SV' (Structural Variant)"
      case symbol_kind::S_ON_INDEL: // "status 'InDel' (Insertion-Deletion)"
        value.copy< kim::VariantFilter::FilterOperation > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_num_op: // num_op
        value.copy< kim::VariantFilter::NumericalOperator > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_str_op: // str_op
        value.copy< kim::VariantFilter::StringOperator > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_STR_VALUE: // "string value"
      case symbol_kind::S_INFO_FIELD: // "info key"
      case symbol_kind::S_string: // string
        value.copy< std::string > (YY_MOVE (that.value));
        break;

      default:
        break;
    }

  }




  template <typename Base>
  VariantFilterParser::symbol_kind_type
  VariantFilterParser::basic_symbol<Base>::type_get () const YY_NOEXCEPT
  {
    return this->kind ();
  }


  template <typename Base>
  bool
  VariantFilterParser::basic_symbol<Base>::empty () const YY_NOEXCEPT
  {
    return this->kind () == symbol_kind::S_YYEMPTY;
  }

  template <typename Base>
  void
  VariantFilterParser::basic_symbol<Base>::move (basic_symbol& s)
  {
    super_type::move (s);
    switch (this->kind ())
    {
      case symbol_kind::S_BOOL_VALUE: // "boolean value"
      case symbol_kind::S_filter: // filter
      case symbol_kind::S_simple_expr: // simple_expr
      case symbol_kind::S_boolean: // boolean
        value.move< bool > (YY_MOVE (s.value));
        break;

      case symbol_kind::S_NUM_VALUE: // "number value"
      case symbol_kind::S_number: // number
        value.move< double > (YY_MOVE (s.value));
        break;

      case symbol_kind::S_INT_VALUE: // "integer value"
      case symbol_kind::S_integer: // integer
        value.move< int > (YY_MOVE (s.value));
        break;

      case symbol_kind::S_bool_op: // bool_op
        value.move< kim::VariantFilter::BooleanOperator > (YY_MOVE (s.value));
        break;

      case symbol_kind::S_EQUAL: // "equality operator ('=')"
      case symbol_kind::S_REGEX: // "regular expression match operator ('~')"
      case symbol_kind::S_NOT_EQUAL: // "inequality operator ('!=' or '<>')"
      case symbol_kind::S_LESS_THAN: // "order relationship operator '<'"
      case symbol_kind::S_LESS_OR_EQUAL: // "order relationship operator '<='"
      case symbol_kind::S_GREATER_THAN: // "order relationship operator '>'"
      case symbol_kind::S_GREATER_OR_EQUAL: // "order relationship operator '>='"
      case symbol_kind::S_any_op: // any_op
        value.move< kim::VariantFilter::ComparisonOperator > (YY_MOVE (s.value));
        break;

      case symbol_kind::S_AND: // "conjunction operator ('AND' or '&&')"
      case symbol_kind::S_OR: // "disjunction operator ('or' or '||')"
      case symbol_kind::S_NOT: // "negation operator ('not' or '!')"
      case symbol_kind::S_BLOCK_OPEN: // "opening parenthesis '('"
      case symbol_kind::S_BLOCK_CLOSE: // "closing parenthesis ')'"
        value.move< kim::VariantFilter::ExpressionOperator > (YY_MOVE (s.value));
        break;

      case symbol_kind::S_ON_CHROM: // "field attribute 'CHROM'"
      case symbol_kind::S_ON_POS: // "field attribute 'POS'"
      case symbol_kind::S_ON_ID: // "field attribute 'ID'"
      case symbol_kind::S_ON_QUAL: // "field attribute 'QUAL'"
      case symbol_kind::S_ON_FILTER: // "field attribute 'FILTER'"
      case symbol_kind::S_ON_INFO: // "field attribute 'INFO'"
      case symbol_kind::S_ON_PLOIDY: // "status 'Ploidy'"
      case symbol_kind::S_ON_SNP: // "status 'SNP' (Single Nucleotide Polymorphism)"
      case symbol_kind::S_ON_MULTI_ALLELIC_SNP: // "status 'MSNP' (multi-allelic SNP)"
      case symbol_kind::S_ON_SV: // "status 'SV' (Structural Variant)"
      case symbol_kind::S_ON_INDEL: // "status 'InDel' (Insertion-Deletion)"
        value.move< kim::VariantFilter::FilterOperation > (YY_MOVE (s.value));
        break;

      case symbol_kind::S_num_op: // num_op
        value.move< kim::VariantFilter::NumericalOperator > (YY_MOVE (s.value));
        break;

      case symbol_kind::S_str_op: // str_op
        value.move< kim::VariantFilter::StringOperator > (YY_MOVE (s.value));
        break;

      case symbol_kind::S_STR_VALUE: // "string value"
      case symbol_kind::S_INFO_FIELD: // "info key"
      case symbol_kind::S_string: // string
        value.move< std::string > (YY_MOVE (s.value));
        break;

      default:
        break;
    }

    location = YY_MOVE (s.location);
  }

  // by_kind.
  inline
  VariantFilterParser::by_kind::by_kind () YY_NOEXCEPT
    : kind_ (symbol_kind::S_YYEMPTY)
  {}

#if 201103L <= YY_CPLUSPLUS
  inline
  VariantFilterParser::by_kind::by_kind (by_kind&& that) YY_NOEXCEPT
    : kind_ (that.kind_)
  {
    that.clear ();
  }
#endif

  inline
  VariantFilterParser::by_kind::by_kind (const by_kind& that) YY_NOEXCEPT
    : kind_ (that.kind_)
  {}

  inline
  VariantFilterParser::by_kind::by_kind (token_kind_type t) YY_NOEXCEPT
    : kind_ (yytranslate_ (t))
  {}



  inline
  void
  VariantFilterParser::by_kind::clear () YY_NOEXCEPT
  {
    kind_ = symbol_kind::S_YYEMPTY;
  }

  inline
  void
  VariantFilterParser::by_kind::move (by_kind& that)
  {
    kind_ = that.kind_;
    that.clear ();
  }

  inline
  VariantFilterParser::symbol_kind_type
  VariantFilterParser::by_kind::kind () const YY_NOEXCEPT
  {
    return kind_;
  }


  inline
  VariantFilterParser::symbol_kind_type
  VariantFilterParser::by_kind::type_get () const YY_NOEXCEPT
  {
    return this->kind ();
  }


#line 15 "variant_filter_parser.ypp"
} // kim
#line 2522 "variant_filter_parser.h"




#endif // !YY_YY_VARIANT_FILTER_PARSER_H_INCLUDED
