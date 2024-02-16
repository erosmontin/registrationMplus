	

ParametersType  readInitFromTxt(RegistrationType::Pointer r, ParametersType initialParameters, int TIN, char * argv[]  ){

#if (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 5) || ITK_VERSION_MAJOR > 4
  itk::TransformFileReaderTemplate<double>::Pointer reader =
    itk::TransformFileReaderTemplate<double>::New();
#else
  itk::TransformFileReader::Pointer reader = itk::TransformFileReader::New();
#endif


  reader->SetFileName(argv[TIN]);
  reader->Update();
  initialParameters= reader->GetTransformList()->begin()->GetPointer()->GetParameters();
  return initialParameters;
};


void writeTranformationTxt(TransformType::Pointer t, int TOUT, char * argv[]  ){
#if (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 5) || ITK_VERSION_MAJOR > 4
				itk::TransformFileWriterTemplate<double>::Pointer writer =
						itk::TransformFileWriterTemplate<double>::New();
	#else
				itk::TransformFileWriter::Pointer writer = itk::TransformFileWriter::New();
	#endif

				//writer->SetInput(registration->GetOutput()->Get());
				writer->SetInput(t);
				writer->SetFileName(argv[TOUT]);
				writer->Update();

};
