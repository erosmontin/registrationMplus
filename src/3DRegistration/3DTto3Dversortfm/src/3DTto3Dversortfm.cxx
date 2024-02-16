#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"

#include "itkTranslationTransform.h"
#include "itkVersorRigid3DTransform.h"

#include "itkTransformFactory.h"


int main( int argc, char * argv[] )
{
  if( argc < 2 )
    {
    std::cerr << "Usage: " << argv[0] << " transformFile" << std::endl;
    return EXIT_FAILURE;
    }
  const char * transformFileName = argv[1];
  typedef double ScalarType;
  const unsigned int Dimension = 3;






  typedef itk::TranslationTransform< ScalarType, Dimension > TranslationTransformType;
  TranslationTransformType::Pointer translation = TranslationTransformType::New();


TranslationTransformType::ParametersType translationParameters;
#if (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 5) || ITK_VERSION_MAJOR > 4
  itk::TransformFileReaderTemplate<double>::Pointer readerT =
    itk::TransformFileReaderTemplate<double>::New();
#else
  itk::TransformFileReader::Pointer readerT = itk::TransformFileReader::New();
#endif

  int TIN=1;
  readerT->SetFileName(argv[TIN]);
  readerT->Update();
  translationParameters= readerT->GetTransformList()->begin()->GetPointer()->GetParameters();
  translation->SetParameters(translationParameters);

  std::cout<<translationParameters<<std::endl;

  typedef itk::VersorRigid3DTransform< ScalarType>VersorRigidTransformType;
  VersorRigidTransformType::Pointer rigid = VersorRigidTransformType::New();


  VersorRigidTransformType::ParametersType rigidParameters;
#if (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 5) || ITK_VERSION_MAJOR > 4
  itk::TransformFileReaderTemplate<double>::Pointer readerR =
    itk::TransformFileReaderTemplate<double>::New();
#else
  itk::TransformFileReader::Pointer readerR = itk::TransformFileReader::New();
#endif

  int RIN=2;
  readerR->SetFileName(argv[RIN]);
  readerR->Update();
  rigidParameters= readerR->GetTransformList()->begin()->GetPointer()->GetParameters();


  rigid->SetParameters(rigidParameters);

  std::cout<<rigidParameters<<std::endl;


  int TOUT=3;

VersorRigidTransformType::ParametersType W;

W=rigidParameters;
for (int p=0;p<3;p++)
{W[p+3]=translationParameters[p];};




VersorRigidTransformType::Pointer rt=VersorRigidTransformType::New();

rt->SetParameters(W);

rt->SetFixedParameters(rigid->GetFixedParameters());

		std::cout<<W<<std::endl;

#if (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 5) || ITK_VERSION_MAJOR > 4
				itk::TransformFileWriterTemplate<double>::Pointer writer =
						itk::TransformFileWriterTemplate<double>::New();
	#else
				itk::TransformFileWriter::Pointer writer = itk::TransformFileWriter::New();
	#endif

				writer->SetInput(rt);
				writer->SetFileName(argv[TOUT]);
				writer->Update();

				std::cout<<rt->GetFixedParameters()<<std::endl;
				std::cout<<rt->GetMatrix()<<std::endl;
				std::cout<<rt->GetTranslation()<<std::endl;
				std::cout<<rt->GetCenter()<<std::endl;


  return EXIT_SUCCESS;
}
