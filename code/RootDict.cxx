// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME codedIRootDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "code/Waveform.h"
#include "code/LecroyFile.h"
#include "code/V1730File.h"
#include "code/TextFile.h"
#include "code/WaveformProcessor.h"
#include "code/DataFile.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static TClass *Waveform_Dictionary();
   static void Waveform_TClassManip(TClass*);
   static void delete_Waveform(void *p);
   static void deleteArray_Waveform(void *p);
   static void destruct_Waveform(void *p);
   static void directoryAutoAdd_Waveform(void *obj, TDirectory *dir);
   static Long64_t merge_Waveform(void *obj, TCollection *coll,TFileMergeInfo *info);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Waveform*)
   {
      ::Waveform *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Waveform));
      static ::ROOT::TGenericClassInfo 
         instance("Waveform", "code/Waveform.h", 8,
                  typeid(::Waveform), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Waveform_Dictionary, isa_proxy, 0,
                  sizeof(::Waveform) );
      instance.SetDelete(&delete_Waveform);
      instance.SetDeleteArray(&deleteArray_Waveform);
      instance.SetDestructor(&destruct_Waveform);
      instance.SetDirectoryAutoAdd(&directoryAutoAdd_Waveform);
      instance.SetMerge(&merge_Waveform);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Waveform*)
   {
      return GenerateInitInstanceLocal((::Waveform*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::Waveform*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Waveform_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::Waveform*)0x0)->GetClass();
      Waveform_TClassManip(theClass);
   return theClass;
   }

   static void Waveform_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *LecroyFile_Dictionary();
   static void LecroyFile_TClassManip(TClass*);
   static void delete_LecroyFile(void *p);
   static void deleteArray_LecroyFile(void *p);
   static void destruct_LecroyFile(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::LecroyFile*)
   {
      ::LecroyFile *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::LecroyFile));
      static ::ROOT::TGenericClassInfo 
         instance("LecroyFile", "code/LecroyFile.h", 25,
                  typeid(::LecroyFile), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &LecroyFile_Dictionary, isa_proxy, 0,
                  sizeof(::LecroyFile) );
      instance.SetDelete(&delete_LecroyFile);
      instance.SetDeleteArray(&deleteArray_LecroyFile);
      instance.SetDestructor(&destruct_LecroyFile);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::LecroyFile*)
   {
      return GenerateInitInstanceLocal((::LecroyFile*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::LecroyFile*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *LecroyFile_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::LecroyFile*)0x0)->GetClass();
      LecroyFile_TClassManip(theClass);
   return theClass;
   }

   static void LecroyFile_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *V1730File_Dictionary();
   static void V1730File_TClassManip(TClass*);
   static void delete_V1730File(void *p);
   static void deleteArray_V1730File(void *p);
   static void destruct_V1730File(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::V1730File*)
   {
      ::V1730File *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::V1730File));
      static ::ROOT::TGenericClassInfo 
         instance("V1730File", "code/V1730File.h", 24,
                  typeid(::V1730File), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &V1730File_Dictionary, isa_proxy, 0,
                  sizeof(::V1730File) );
      instance.SetDelete(&delete_V1730File);
      instance.SetDeleteArray(&deleteArray_V1730File);
      instance.SetDestructor(&destruct_V1730File);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::V1730File*)
   {
      return GenerateInitInstanceLocal((::V1730File*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::V1730File*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *V1730File_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::V1730File*)0x0)->GetClass();
      V1730File_TClassManip(theClass);
   return theClass;
   }

   static void V1730File_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *WaveformProcessor_Dictionary();
   static void WaveformProcessor_TClassManip(TClass*);
   static void delete_WaveformProcessor(void *p);
   static void deleteArray_WaveformProcessor(void *p);
   static void destruct_WaveformProcessor(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WaveformProcessor*)
   {
      ::WaveformProcessor *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::WaveformProcessor));
      static ::ROOT::TGenericClassInfo 
         instance("WaveformProcessor", "code/WaveformProcessor.h", 35,
                  typeid(::WaveformProcessor), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &WaveformProcessor_Dictionary, isa_proxy, 0,
                  sizeof(::WaveformProcessor) );
      instance.SetDelete(&delete_WaveformProcessor);
      instance.SetDeleteArray(&deleteArray_WaveformProcessor);
      instance.SetDestructor(&destruct_WaveformProcessor);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WaveformProcessor*)
   {
      return GenerateInitInstanceLocal((::WaveformProcessor*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::WaveformProcessor*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *WaveformProcessor_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::WaveformProcessor*)0x0)->GetClass();
      WaveformProcessor_TClassManip(theClass);
   return theClass;
   }

   static void WaveformProcessor_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrapper around operator delete
   static void delete_Waveform(void *p) {
      delete ((::Waveform*)p);
   }
   static void deleteArray_Waveform(void *p) {
      delete [] ((::Waveform*)p);
   }
   static void destruct_Waveform(void *p) {
      typedef ::Waveform current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around the directory auto add.
   static void directoryAutoAdd_Waveform(void *p, TDirectory *dir) {
      ((::Waveform*)p)->DirectoryAutoAdd(dir);
   }
   // Wrapper around the merge function.
   static Long64_t  merge_Waveform(void *obj,TCollection *coll,TFileMergeInfo *) {
      return ((::Waveform*)obj)->Merge(coll);
   }
} // end of namespace ROOT for class ::Waveform

namespace ROOT {
   // Wrapper around operator delete
   static void delete_LecroyFile(void *p) {
      delete ((::LecroyFile*)p);
   }
   static void deleteArray_LecroyFile(void *p) {
      delete [] ((::LecroyFile*)p);
   }
   static void destruct_LecroyFile(void *p) {
      typedef ::LecroyFile current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::LecroyFile

namespace ROOT {
   // Wrapper around operator delete
   static void delete_V1730File(void *p) {
      delete ((::V1730File*)p);
   }
   static void deleteArray_V1730File(void *p) {
      delete [] ((::V1730File*)p);
   }
   static void destruct_V1730File(void *p) {
      typedef ::V1730File current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::V1730File

namespace ROOT {
   // Wrapper around operator delete
   static void delete_WaveformProcessor(void *p) {
      delete ((::WaveformProcessor*)p);
   }
   static void deleteArray_WaveformProcessor(void *p) {
      delete [] ((::WaveformProcessor*)p);
   }
   static void destruct_WaveformProcessor(void *p) {
      typedef ::WaveformProcessor current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WaveformProcessor

namespace {
  void TriggerDictionaryInitialization_RootDict_Impl() {
    static const char* headers[] = {
"code/Waveform.h",
"code/LecroyFile.h",
"code/V1730File.h",
"code/TextFile.h",
"code/WaveformProcessor.h",
"code/DataFile.h",
0
    };
    static const char* includePaths[] = {
"/include",
"/Users/luca/miniforge3/include/",
"/Users/luca/Giacomo/WFprocessor/luca/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "RootDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$code/Waveform.h")))  Waveform;
class __attribute__((annotate("$clingAutoload$code/LecroyFile.h")))  LecroyFile;
class __attribute__((annotate("$clingAutoload$code/V1730File.h")))  V1730File;
class __attribute__((annotate("$clingAutoload$code/WaveformProcessor.h")))  WaveformProcessor;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "RootDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "code/Waveform.h"
#include "code/LecroyFile.h"
#include "code/V1730File.h"
#include "code/TextFile.h"
#include "code/WaveformProcessor.h"
#include "code/DataFile.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"LecroyFile", payloadCode, "@",
"V1730File", payloadCode, "@",
"Waveform", payloadCode, "@",
"WaveformProcessor", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("RootDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_RootDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_RootDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_RootDict() {
  TriggerDictionaryInitialization_RootDict_Impl();
}
