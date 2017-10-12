TGraphErrors* ReadHepdata(const char* fileName, Float_t normalizeTo = 0, Bool_t errorsAreAdded = kFALSE, Int_t skipYErrors = 0, Int_t skipXerrors = 1)
{
  // expected format: x [x2] y [ye] [xe]
  //
  // skipYErrors:   0 --> ye present
  //                1 --> no errors ye
  //                2 --> y and ye are lower and upper error, i.e. y' = (y + ye) / 2 and ye = (ye - y) / 2
  // 
  // skipXerrors:   0 --> xe present
  //                1 --> no errors xe
  //                2 --> x2 present, xe not present and is calculated from x2 - x
  
  ifstream fin(fileName);

  graph = new TGraphErrors(0);

  Double_t sum = 0;

  while (fin.good())
  {
    char buffer[2000];
    if (fin.peek() == '#')
    {
      fin.getline(buffer, 2000);
      continue;
    }
  
    Float_t x = -1;
    Float_t x2 = -1;
    Float_t y = -1;
    Float_t ye = 0;
    Float_t xe = 0;

    fin >> x;
    
    if (skipXerrors == 2)
    {
      fin >> x2;
      xe = (x2 - x + 1) / 2;
      x = x + (x2 - x) / 2;
    }
    
    fin >> y;

    if (y == -1)
      continue;

    if (skipYErrors == 0)
    {
      ye = -1;
      fin >> ye;
      if (ye == -1)
        continue;
    }
    else if (skipYErrors == 2)
    {
      ye = -1;
      fin >> ye;
      if (ye == -1)
        continue;
      
      Float_t newy = (y + ye) / 2;
      ye = (ye - y) / 2;
      y = newy;
    }

    if (skipXerrors == 0)
    {
      xe = -1;
      fin >> xe;
      if (xe == -1)
        continue;
    }

    //Printf("%f %f %f %f", x, y, xe, ye);

    if (errorsAreAdded)
      ye -= y;

    graph->SetPoint(graph->GetN(), x, y);
    graph->SetPointError(graph->GetN()-1, xe, ye);

    sum += y;
    
    // read rest until end of line...
    fin.getline(buffer, 2000);
  }
  fin.close();

  Printf("%s: %f", fileName, sum);

  NormalizeTo(graph, normalizeTo);

  return graph;
}

TGraphErrors* ReadHepdataInput(const char* fileName, Float_t normalizeTo = 0)
{
  // expected format: xleft TO xright; y +- ye;
  
  ifstream fin(fileName);

  graph = new TGraphErrors(0);

  Double_t sum = 0;

  while (fin.good())
  {
    char buffer[2000];
    fin.getline(buffer, 2000);
    if (buffer[0] == '*')
      continue;

    Float_t x = -1;
    Float_t x2 = -1;
    Float_t y = -1;
    Float_t ye = 0;
    
    sscanf(buffer, "%f TO %f; %f +- %f", &x, &x2, &y, &ye);

    Printf("%f %f %f %f", x, x2, y, ye);

    graph->SetPoint(graph->GetN(), (x + x2) / 2, y);
    graph->SetPointError(graph->GetN()-1, (x2 - x) / 2, ye);

    sum += y;
  }
  fin.close();

  Printf("%s: %f", fileName, sum);

  NormalizeTo(graph, normalizeTo);

  return graph;
}

void ScaleGraph(TGraphErrors* graph, Float_t factor)
{
        for (Int_t i=0; i<graph->GetN(); i++)
        {
    graph->GetY()[i] *= factor;
    graph->GetEY()[i] *= factor;
  }
}

TGraphErrors* CMS()
{
  // Table 14: Fully corrected charged hadron multiplicity spectrum for |pseudorapidity| < 1.5 at a centre-of-mass energy of 7000 GeV.. 
  // http://hepdata.cedar.ac.uk/view/ins879315/d14
  // WARNING errors symmetrized below!

  // Plot: p8068_d14x1y1
  double p8068_d14x1y1_xval[] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 
    9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 
    19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 
    29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 
    39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 
    49.0, 50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 
    59.0, 60.0, 61.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 
    69.0, 70.0, 71.0, 72.0, 73.0, 74.0, 75.0, 76.0, 77.0, 78.0, 
    79.5, 81.5, 83.5, 85.5, 88.0, 91.5, 96.5, 102.5, 108.5, 114.5, 
    120.5, 126.5, 132.5, 138.5, 144.5, 150.5 };
  double p8068_d14x1y1_xerrminus[] = { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
    1.0, 1.0, 1.0, 1.0, 1.5, 2.0, 3.0, 3.0, 3.0, 3.0, 
    3.0, 3.0, 3.0, 3.0, 3.0, 3.0 };
  double p8068_d14x1y1_xerrplus[] = { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
    1.0, 1.0, 1.0, 1.0, 1.5, 2.0, 3.0, 3.0, 3.0, 3.0, 
    3.0, 3.0, 3.0, 3.0, 3.0, 3.0 };
  double p8068_d14x1y1_yval[] = { 0.0604, 0.02425, 0.02784, 0.0317, 0.03515, 0.03769, 0.03928, 0.0396, 0.03881, 
    0.03716, 0.03495, 0.03257, 0.03026, 0.02815, 0.02634, 0.02476, 0.02337, 0.02211, 0.02098, 
    0.01992, 0.01894, 0.01802, 0.01715, 0.01634, 0.01559, 0.01487, 0.01419, 0.01353, 0.01291, 
    0.01231, 0.01174, 0.01119, 0.01064, 0.01009, 0.00956, 0.00905, 0.00858, 0.00815, 0.00774, 
    0.00734, 0.00695, 0.00657, 0.00622, 0.0059, 0.00562, 0.00535, 0.00514, 0.00485, 0.00451, 
    0.00424, 0.00398, 0.00374, 0.00351, 0.00331, 0.00313, 0.00296, 0.00279, 0.00263, 0.00247, 
    0.0023, 0.00215, 0.002005, 0.00187, 0.001743, 0.001627, 0.001525, 0.001431, 0.00134, 0.00125, 
    0.001162, 0.00108, 0.001006, 9.41E-4, 8.83E-4, 8.27E-4, 7.73E-4, 7.19E-4, 6.66E-4, 6.14E-4, 
    5.37E-4, 4.62E-4, 4.03E-4, 3.43E-4, 2.64E-4, 1.96E-4, 1.3E-4, 8.3E-5, 4.83E-5, 2.62E-5, 
    1.39E-5, 7.5E-6, 4.2E-6, 2.27E-6, 1.11E-6, 4.8E-7 };
  double p8068_d14x1y1_yerrminus[] = { 0.01865261375786246, 0.008312412405553516, 0.008805112151472007, 0.008666285248017169, 0.007655494758668443, 0.006125952987086988, 0.004538039224158381, 0.003151586901863885, 0.0027025358462007495, 
    0.0024030397416605494, 0.002064000968992021, 0.0018060177186284745, 0.0016080111939908877, 0.0014500000000000001, 0.0013497407158413798, 0.0012495199078045936, 0.00117085438889727, 0.0010825894882179486, 0.0010022474744293447, 
    9.144397191723465E-4, 8.464632301523795E-4, 7.692203845452875E-4, 7.021395872616784E-4, 6.419501538281613E-4, 5.818934610390462E-4, 5.255473337388365E-4, 4.884669896727925E-4, 4.5188494110780014E-4, 4.1593268686170845E-4, 
    3.9357337308308855E-4, 3.7589892258425004E-4, 3.584689665786984E-4, 3.49857113690718E-4, 3.4132096331752025E-4, 3.27566787083184E-4, 3.106444913401813E-4, 2.968164415931166E-4, 2.8844410203711916E-4, 2.8844410203711916E-4, 
    2.80178514522438E-4, 2.80178514522438E-4, 2.7459060435491964E-4, 2.66270539113887E-4, 2.66270539113887E-4, 2.66270539113887E-4, 2.6076809620810597E-4, 2.8653097563788803E-4, 2.778488797889961E-4, 2.5553864678361276E-4, 
    2.5553864678361276E-4, 2.4698178070456936E-4, 2.4698178070456936E-4, 2.3323807579381201E-4, 2.3323807579381201E-4, 2.3323807579381201E-4, 2.1954498400100151E-4, 2.1954498400100151E-4, 2.1954498400100151E-4, 2.147091055358389E-4, 
    2.0591260281974002E-4, 1.972308292331602E-4, 1.9449421585229727E-4, 1.878004259846074E-4, 1.7887705274852894E-4, 1.6995587662684689E-4, 1.6240997506310995E-4, 1.5657905351610732E-4, 1.5212166183683373E-4, 1.471461858153313E-4, 
    1.4183441049336372E-4, 1.3567977004697496E-4, 1.290038759107648E-4, 1.2317467272130256E-4, 1.1819052415485769E-4, 1.1320777358467925E-4, 1.0960383204979651E-4, 1.0462313319720452E-4, 9.964436762808022E-5, 9.364827814754525E-5, 
    8.452810183601664E-5, 7.026378868236468E-5, 6.276941930590087E-5, 5.7939623747483895E-5, 4.5541190146942804E-5, 3.54400902933387E-5, 2.4478766308782802E-5, 1.7488853593074647E-5, 1.158015543937127E-5, 7.15122367151245E-6, 
    4.248529157249601E-6, 2.523885892824792E-6, 1.6278820596099706E-6, 1.0780074211247343E-6, 7.247068372797376E-7, 4.588027898781785E-7 };
  double p8068_d14x1y1_yerrplus[] = { 0.04482186966202994, 0.008811140675304191, 0.009194895322949577, 0.009016041259887846, 0.008015247968715628, 0.0065155966726002925, 0.005007284693324317, 0.0037497333238511775, 0.003310226578347772, 
    0.0030701954335188497, 0.002810266891240047, 0.0025612692166189794, 0.0023423279018958895, 0.0021434784813475502, 0.0019734487578855447, 0.0018233211456021672, 0.001684428686528462, 0.001565535052306399, 0.0014354441821262157, 
    0.0013267252918370102, 0.001218236430254817, 0.001119866063420086, 0.001031600697944704, 9.414881836751856E-4, 8.612200647918045E-4, 7.93095202355934E-4, 7.253275122315435E-4, 6.676076692189808E-4, 6.198386886924694E-4, 
    5.692099788303082E-4, 5.314132102234569E-4, 4.939635614091387E-4, 5.692099788303082E-4, 4.569463863518345E-4, 4.0718546143004664E-4, 3.8910152916687437E-4, 3.577708763999664E-4, 3.3999999999999997E-4, 3.312099032335839E-4, 
    3.138470965295043E-4, 3.0528675044947496E-4, 2.91547594742265E-4, 2.830194339616981E-4, 2.830194339616981E-4, 2.7459060435491964E-4, 2.6925824035672517E-4, 2.6076809620810597E-4, 2.6076809620810597E-4, 2.729468812791236E-4, 
    2.729468812791236E-4, 2.641968962724581E-4, 2.641968962724581E-4, 2.505992817228334E-4, 2.418677324489565E-4, 2.418677324489565E-4, 2.2825424421026657E-4, 2.2825424421026657E-4, 2.1954498400100151E-4, 2.147091055358389E-4, 
    2.0591260281974002E-4, 2.0591260281974002E-4, 1.9535864454894234E-4, 1.878004259846074E-4, 1.7887705274852894E-4, 1.8124568960391856E-4, 1.7104677722775134E-4, 1.5572411502397438E-4, 1.512679741386127E-4, 1.462908062729849E-4, 
    1.4183441049336372E-4, 1.3652838532700809E-4, 1.306981254647518E-4, 1.2486793023030372E-4, 1.1988744721612851E-4, 1.1576268828944842E-4, 1.1216059914247962E-4, 1.0718675291284833E-4, 1.030776406404415E-4, 9.712878049270463E-5, 
    8.814193099768125E-5, 7.392563831310488E-5, 6.74166151627327E-5, 6.168468205316454E-5, 5.508175741568164E-5, 3.83275357934736E-5, 2.757698315624826E-5, 1.942832983043061E-5, 1.283160161476345E-5, 8.000624975587846E-6, 
    4.884669896727926E-6, 3.04138126514911E-6, 2.0248456731316586E-6, 1.3805795884337853E-6, 9.265527507918802E-7, 5.73846669416143E-7 };
  double p8068_d14x1y1_ystatminus[] = { 0.0014, 6.1E-4, 3.0E-4, 3.3E-4, 2.9E-4, 2.7E-4, 2.7E-4, 2.7E-4, 2.6E-4, 
    2.5E-4, 2.4E-4, 2.4E-4, 2.4E-4, 2.4E-4, 2.3E-4, 2.2E-4, 2.2E-4, 2.2E-4, 2.1E-4, 
    2.1E-4, 2.1E-4, 2.1E-4, 2.1E-4, 2.0E-4, 1.9E-4, 1.9E-4, 1.9E-4, 1.9E-4, 1.9E-4, 
    1.8E-4, 1.8E-4, 1.8E-4, 1.8E-4, 1.8E-4, 1.7E-4, 1.7E-4, 1.6E-4, 1.6E-4, 1.6E-4, 
    1.6E-4, 1.6E-4, 1.5E-4, 1.5E-4, 1.5E-4, 1.5E-4, 1.4E-4, 1.4E-4, 1.4E-4, 1.3E-4, 
    1.3E-4, 1.3E-4, 1.3E-4, 1.2E-4, 1.2E-4, 1.2E-4, 1.1E-4, 1.1E-4, 1.1E-4, 1.0E-4, 
    1.0E-4, 1.0E-4, 9.8E-5, 9.5E-5, 9.1E-5, 8.7E-5, 8.4E-5, 8.1E-5, 7.9E-5, 7.6E-5, 
    7.4E-5, 7.2E-5, 6.9E-5, 6.6E-5, 6.3E-5, 6.0E-5, 5.8E-5, 5.5E-5, 5.2E-5, 4.7E-5, 
    3.7E-5, 2.9E-5, 2.4E-5, 2.1E-5, 1.5E-5, 1.0E-5, 6.5E-6, 4.5E-6, 3.3E-6, 2.5E-6, 
    1.9E-6, 1.4E-6, 1.1E-6, 8.6E-7, 6.4E-7, 4.3E-7 };
  double p8068_d14x1y1_ystatplus[] = { 0.0014, 6.1E-4, 3.0E-4, 3.3E-4, 2.9E-4, 2.7E-4, 2.7E-4, 2.7E-4, 2.6E-4, 
    2.5E-4, 2.4E-4, 2.4E-4, 2.4E-4, 2.4E-4, 2.3E-4, 2.2E-4, 2.2E-4, 2.2E-4, 2.1E-4, 
    2.1E-4, 2.1E-4, 2.1E-4, 2.1E-4, 2.0E-4, 1.9E-4, 1.9E-4, 1.9E-4, 1.9E-4, 1.9E-4, 
    1.8E-4, 1.8E-4, 1.8E-4, 1.8E-4, 1.8E-4, 1.7E-4, 1.7E-4, 1.6E-4, 1.6E-4, 1.6E-4, 
    1.6E-4, 1.6E-4, 1.5E-4, 1.5E-4, 1.5E-4, 1.5E-4, 1.4E-4, 1.4E-4, 1.4E-4, 1.3E-4, 
    1.3E-4, 1.3E-4, 1.3E-4, 1.2E-4, 1.2E-4, 1.2E-4, 1.1E-4, 1.1E-4, 1.1E-4, 1.0E-4, 
    1.0E-4, 1.0E-4, 9.8E-5, 9.5E-5, 9.1E-5, 8.7E-5, 8.4E-5, 8.1E-5, 7.9E-5, 7.6E-5, 
    7.4E-5, 7.2E-5, 6.9E-5, 6.6E-5, 6.3E-5, 6.0E-5, 5.8E-5, 5.5E-5, 5.2E-5, 4.7E-5, 
    3.7E-5, 2.9E-5, 2.4E-5, 2.1E-5, 1.5E-5, 1.0E-5, 6.5E-6, 4.5E-6, 3.3E-6, 2.5E-6, 
    1.9E-6, 1.4E-6, 1.1E-6, 8.6E-7, 6.4E-7, 4.3E-7 };
  int p8068_d14x1y1_numpoints = 95;
  p8068_d14x1y1 = new TGraphErrors(p8068_d14x1y1_numpoints, p8068_d14x1y1_xval, p8068_d14x1y1_yval, p8068_d14x1y1_xerrminus, p8068_d14x1y1_yerrminus);

  return p8068_d14x1y1;
}

void Scale(TGraphErrors* graph, Float_t factor)
{
  for (Int_t i=0; i<graph->GetN(); i++)
  {
    graph->SetPoint(i, graph->GetX()[i],  graph->GetY()[i] * factor);
    graph->SetPointError(i, graph->GetEX()[i],  graph->GetEY()[i] * factor);
  }
}

void NormalizeTo(TGraphErrors* graph, Float_t normalizeTo)
{
        Float_t sum = 0;
        for (Int_t i=0; i<graph->GetN(); i++)
                sum += graph->GetY()[i];
        
        if (normalizeTo > 0 && sum > 0)
        {
                Scale(graph, normalizeTo / sum);
        }       
}

void NormalizeToWidth(TGraphErrors* graph, Float_t normalizeTo)
{
        Float_t sum = 0;
        for (Int_t i=0; i<graph->GetN(); i++)
        {
                Float_t width = 0;
                if (i > 0)
                        width += (graph->GetX()[i] - graph->GetX()[i-1]) / 2;
                if (i < graph->GetN()-1)
                        width += (graph->GetX()[i+1] - graph->GetX()[i]) / 2;
                if (i == 0 || i == graph->GetN() - 1)
                        width *= 2;
                                        
                sum += graph->GetY()[i] * width;
        }
        
        if (normalizeTo > 0 && sum > 0)
        {
                sum /= normalizeTo;
                for (Int_t i=0; i<graph->GetN(); i++)
                {
                        graph->SetPoint(i, graph->GetX()[i],  graph->GetY()[i] / sum);
                        graph->SetPointError(i, graph->GetEX()[i],  graph->GetEY()[i] / sum);
                }
        }       
}

void NormalizeToXError(TGraphErrors* graph, Float_t normalizeTo)
{
        Float_t sum = 0;
        for (Int_t i=0; i<graph->GetN(); i++)
        {
                Float_t width = graph->GetEX()[i] * 2;
                sum += graph->GetY()[i] * width;
        }
        
        if (normalizeTo > 0 && sum > 0)
        {
                sum /= normalizeTo;
                for (Int_t i=0; i<graph->GetN(); i++)
                {
                        graph->SetPoint(i, graph->GetX()[i],  graph->GetY()[i] / sum);
                        graph->SetPointError(i, graph->GetEX()[i],  graph->GetEY()[i] / sum);
                }
        }       
}

TF1* GetNBDLog(const char* name = "nbd")
{
        TF1* func = new TF1(name, "exp(log([0]) + TMath::LnGamma([2]+x) - TMath::LnGamma([2]) - TMath::LnGamma(x+1) + log([1] / ([1]+[2])) * x + log(1.0 + [1]/[2]) * -[2])");
        
        func->SetParNames("scaling", "averagen", "k");
        func->SetParLimits(0, 0.5, 2);
        func->SetParLimits(1, 1, 100);
        func->SetParLimits(2, 1, 20);
        func->SetParameters(1, 10, 2);
        
        return func;
}

void LR()
{
//   k = new TF1("k", "1.0/([0]+[1]*log(x))", 1, 1e5);
//   k->SetParameters(-1.89842e+01/200, 1.91776e+01/200); // fit whole range

  //   k->SetParameters(4.16930e+00/200, 1.66523e+01/200);  // fit above 10

  // fit above 30% of x range
  k = new TF1("k", "1.94", 1, 1e5); 
  n = new TF1("n", "[0]*x**(2*[1])", 25, 100000);
  n->SetParameters(1.636, 0.1436);
  
  // fit in full range
//   k2 = new TF1("k2", "[0]+[1]*log(x)", 25, 100000);
//   k2->SetParameters(-2.267e1/200, 1.909e1/200);
//   n2 = new TF1("n2", "[0]*x**(2*[1])", 25, 100000);
//   n2->SetParameters(2.146, 0.1218);

  // fit above 10% of x range
  k2 = new TF1("k2", "[0]+[1]*log(x)", 25, 100000);
  k2->SetParameters(2.48e1/200, 1.291e1/200);
  n2 = new TF1("n2", "[0]*x**(2*[1])", 25, 100000);
  n2->SetParameters(2.045, 0.1245);

  canvas = new TCanvas("c", "c", 800, 600);
	gPad->SetRightMargin(0.05);
	gPad->SetTopMargin(0.05);
  
  dummy = new TH2F("dummy", ";N_{ch};P(N_{ch})", 100, 0, 400, 100, 1e-7, 0.1);
  dummy->SetStats(0);
	dummy->GetYaxis()->SetTitleOffset(1.2);
  dummy->Draw();
  gPad->SetLogy();
  legend = new TLegend(0.58, 0.6, 0.91, 0.92);
  legend->SetFillColor(0);
	legend->SetLineColor(0);
  
  Int_t count = 6;
  const char* files[] = { "mult_eta15_nsd_200.txt", "mult_eta15_nsd_540.txt", 
  "mult_eta15_nsd_900.txt", 
  "CMS", "fig14b.txt", "fig14c.txt" };
  const char* names[] = { "UA5 200 GeV", "UA5 540 GeV", "UA5 900 GeV", "CMS 7 TeV", "ALICE 2.76 TeV", "ALICE 7 TeV" };
  Float_t energy[] = { 200, 540, 900, 7000, 2760, 7000 };
  Int_t fileType[] = { 0, 2, 0, -1, -2, -2 };
  Int_t norm[] = { 1, 1, 1, 1, 1, 1 }; // normalization to width necessary
  Float_t reduceError[] = { 0, 0, 0, 0, 0, 0, 0, 0 };
  Int_t markers[] = { 20, 21, 22, 23, 24, 25, 25, 26, 27, 28, 30 };

   TFile*fout=new TFile("outputALICE7TeV.root","recreate");

  for (Int_t i=0; i<count; i++)
  {
    if (i == 1 || i == 3)
      continue;
    
    if (fileType[i] == -1)
      graph = CMS();
    else if (fileType[i] == -2)
      graph = ReadHepdataInput(files[i], 1);
    else
      graph = ReadHepdata(files[i], 1, kFALSE, kFALSE, fileType[i]);
    if (norm[i] == 1)
      NormalizeToXError(graph, 1);
    if (norm[i] == 2)
      NormalizeToWidth(graph, 1);
      
    graph->SetMarkerStyle(markers[i+4]);
//     graph->SetMarkerSize(0.8);
    
    graph->Clone()->Draw("PSAME");
    

    if (i==5) { 
      graph->SetName("graphALICE7TeV");
      graph->Write();
   }
    legend->AddEntry(graph, names[i], "P");

  }
  fout->Write();

  
  func = GetNBDLog("nbd_0.2");
  func->SetParameters(1, 8, k->Eval(200));
  func->SetRange(0.3 * 38, 60);
  func->Draw("SAME");
// 	func->SetLineStyle(2);
  legend->AddEntry(func, "NBD 0.2 TeV", "L");
  
  Printf("%f", func->Integral(0, 600));
 
  func = GetNBDLog("nbd_7");
  func->SetParameters(1, 20.8, k->Eval(7000));
  func->SetRange(0.3 * 151, 170);
  func->SetLineColor(2);
// 	func->SetLineStyle(3);
  func->Draw("SAME");
  legend->AddEntry(func, "NBD 7 TeV", "L");
  
  Printf("%f", func->Integral(0, 600));

  // TODO change cms energy here
  
  func = GetNBDLog("nbd_100");
  func->SetParameters(1, n->Eval(100000), k->Eval(100000));
//   func->SetParameters(1, 19 * 1.8, TMath::Max(1.0, k->Eval(100000)));
  func->SetRange(1, 340);
  func->SetLineColor(4);
  func->Draw("SAME");
  legend->AddEntry(func, "NBD 100 TeV", "L");
  
  func = GetNBDLog("nbd_100_2");
  func->SetParameters(1, n2->Eval(100000), k2->Eval(100000));
  func->SetRange(1, 800);
  func->SetLineColor(4);
  func->SetLineStyle(2);
  //func->Draw("SAME");
  //legend->AddEntry(func, "NBD 100 TeV", "L");

  Printf("%f", func->Integral(0, 600));

  legend->Draw();  

  canvas->SaveAs("fcc_eta15.png");  
}
