void analysisFJ(){

    // Load fastjet libraries and includes
    gSystem->AddIncludePath( "-I$FASTJET/include");
    gSystem->Load("$FASTJET/lib/libfastjet");

    gROOT->ProcessLine(".L Test.C");
    gROOT->ProcessLine("Test t;");
    gROOT->ProcessLine("t.Loop()");
}