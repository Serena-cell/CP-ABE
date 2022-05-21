import it.unisa.dia.gas.jpbc.Element;
import it.unisa.dia.gas.jpbc.Pairing;
import it.unisa.dia.gas.plaf.jpbc.pairing.PairingFactory;

import java.io.*;
import java.math.BigInteger;
import java.util.*;



public class CPABE {

    //generate random number ∈ Zr alpha, beta
    public static void setup(String pairingParametersFileName, String pkFileName, String mkFileName) {
        Pairing bp = PairingFactory.getPairing(pairingParametersFileName);

        Element g = bp.getG1().newRandomElement().getImmutable();
        Element alpha = bp.getZr().newRandomElement().getImmutable();
        Element beta = bp.getZr().newRandomElement().getImmutable();
        Element h=bp.getGT().newRandomElement().getImmutable();

        Element g_alpha = g.powZn(alpha).getImmutable();
        Element g_beta = g.powZn(beta).getImmutable();
        Element g_alpha_beta=(g.powZn(alpha)).powZn(beta);
        Element egg_alpha =(bp.pairing(g,g).powZn(alpha)).getImmutable();
        Element egg_alpha_beta = bp.pairing(g,g).powZn(alpha.mul(beta)).getImmutable();

        //MK={beta,g_alpha}
        Properties mkProp = new Properties();
        mkProp.setProperty("g_alpha", Base64.getEncoder().withoutPadding().encodeToString(g_alpha.toBytes()));
        mkProp.setProperty("beta", Base64.getEncoder().withoutPadding().encodeToString(beta.toBytes()));

        //PK={g,h,g_beta,g_alpha_beta,egg_alpha,egg_alpha_beta}
        Properties pkProp = new Properties();
        pkProp.setProperty("g", Base64.getEncoder().withoutPadding().encodeToString(g.toBytes()));
        pkProp.setProperty("h", Base64.getEncoder().withoutPadding().encodeToString(h.toBytes()));
       // pkProp.setProperty("g_alpha", Base64.getEncoder().withoutPadding().encodeToString(g_alpha.toBytes()));
        pkProp.setProperty("g_beta", Base64.getEncoder().withoutPadding().encodeToString(g_beta.toBytes()));
        pkProp.setProperty("g_alpha_beta", Base64.getEncoder().withoutPadding().encodeToString(g_alpha_beta.toBytes()));
        pkProp.setProperty("egg_alpha", Base64.getEncoder().withoutPadding().encodeToString(egg_alpha.toBytes()));
        pkProp.setProperty("egg_alpha_beta", Base64.getEncoder().withoutPadding().encodeToString(egg_alpha_beta.toBytes()));

        storePropToFile(mkProp, mkFileName);
        storePropToFile(pkProp, pkFileName);
    }

    // to prepare for later keygen, identical to the attribute authority setup in my paper
    // here we suppose that an authority holds all attributes from the input authorityAttList
    //generate random number ∈ Zr t, delta_att
    public static void aa_keygen(String pairingParametersFileName, int[] authorityAttList, String pkFileName, String mkFileName,
                                 String sk_aaFileName){
        Pairing bp = PairingFactory.getPairing(pairingParametersFileName);
        Properties sk_aaProp = loadPropFromFile(sk_aaFileName);
        Properties pkProp = loadPropFromFile(pkFileName);

        //resolve information from pkFileName
        String gString = pkProp.getProperty("g");
        Element g = bp.getG1().newElementFromBytes(Base64.getDecoder().decode(gString)).getImmutable();

        Properties mkProp = loadPropFromFile(mkFileName);
        //resolve information from mkFileName
        String g_alphaString = mkProp.getProperty("g_alpha");
        Element g_alpha = bp.getG1().newElementFromBytes(Base64.getDecoder().decode(g_alphaString)).getImmutable();

//        List<Integer> accessTreeAttList = new ArrayList<>();
//        int[] accessTreeAtt;
//        int j;
//        for(j=0;j<accessTree.length;j++)
//        {
//            //Extract the leaf node information of the visit tree
//            if(accessTree[j].isLeaf()){
//                accessTreeAttList.add(accessTree[j].att);
//            }
//        }
//        accessTreeAtt = accessTreeAttList.stream().mapToInt(i->i).toArray();
        //t equals gama in my paper, as the algorithm shows what happen in one authority.
        Element t = bp.getZr().newRandomElement().getImmutable();

        for (int att:authorityAttList) {
            //for each attribute in the accessTreeAttList:
            Element delta_att=bp.getZr().newRandomElement().getImmutable();
            Element D_aa =g_alpha.powZn((delta_att.mul(t)).pow(BigInteger.valueOf(-1)));
            Element D_aa_=g.powZn(delta_att).powZn(t);
            sk_aaProp.setProperty("D_aa"+att, Base64.getEncoder().withoutPadding().encodeToString(D_aa.toBytes()));
            sk_aaProp.setProperty("D_aa_"+att, Base64.getEncoder().withoutPadding().encodeToString(D_aa_.toBytes()));
        }
        storePropToFile(sk_aaProp, sk_aaFileName);
    }

    //generate random number ∈ Zr theta, tau
    public static void keygen(String pairingParametersFileName, int[] userAttList, String pkFileName, String mkFileName,
                              String sk_aaFileName, String skFileName) {
        Pairing bp = PairingFactory.getPairing(pairingParametersFileName);

        Properties pkProp = loadPropFromFile(pkFileName);
        Properties sk_aaProp=loadPropFromFile(sk_aaFileName);

        //resolve information from pkFileName
        String gString = pkProp.getProperty("g");
        Element g = bp.getG1().newElementFromBytes(Base64.getDecoder().decode(gString)).getImmutable();
        String hString = pkProp.getProperty("h");
        Element h = bp.getG1().newElementFromBytes(Base64.getDecoder().decode(hString)).getImmutable();
        String g_alpha_betaString = pkProp.getProperty("g_alpha_beta");
        Element g_alpha_beta = bp.getG1().newElementFromBytes(Base64.getDecoder().decode(g_alpha_betaString)).getImmutable();
        String egg_alphaString = pkProp.getProperty("egg_alpha");
        Element egg_alpha = bp.getGT().newElementFromBytes(Base64.getDecoder().decode(egg_alphaString)).getImmutable();

        Properties mkProp = loadPropFromFile(mkFileName);
        //resolve information from mkFileName
        String g_alphaString = mkProp.getProperty("g_alpha");
        Element g_alpha = bp.getG1().newElementFromBytes(Base64.getDecoder().decode(g_alphaString)).getImmutable();

        Properties skProp = new Properties();

        Element theta = bp.getZr().newRandomElement().getImmutable();
        Element tau = bp.getZr().newRandomElement().getImmutable();

        //D3=egg_alpha_theta
        Element D3=egg_alpha.powZn(theta);
        Element D2=g.powZn(tau).getImmutable();
        Element h_tau=h.powZn(tau).getImmutable();
        Element D1 = (g_alpha.powZn(theta)).mul(h_tau).getImmutable();
        Element D=g_alpha_beta.mul(g_alpha.powZn(theta));

        skProp.setProperty("D", Base64.getEncoder().withoutPadding().encodeToString(D.toBytes()));
        skProp.setProperty("D1", Base64.getEncoder().withoutPadding().encodeToString(D1.toBytes()));
        skProp.setProperty("D2", Base64.getEncoder().withoutPadding().encodeToString(D2.toBytes()));
        skProp.setProperty("D3", Base64.getEncoder().withoutPadding().encodeToString(D3.toBytes()));

        //hint: att in userAttList must be a member of authorityAttList
        for (int att:userAttList ) {
            //resolve information from sk_aaFileName
            String D_aaString = sk_aaProp.getProperty("D_aa"+att);
            Element D_aa = bp.getG1().newElementFromBytes(Base64.getDecoder().decode(D_aaString)).getImmutable();
            Element D_usr =D_aa.powZn(theta);
            skProp.setProperty("D_usr"+att, Base64.getEncoder().withoutPadding().encodeToString(D_usr.toBytes()));
        }

        skProp.setProperty("userAttList", Arrays.toString(userAttList));
        storePropToFile(skProp, skFileName);
    }

    //generate random number ∈ Zr omega, s
    public static void encrypt(String pairingParametersFileName, Element message, Node[] accessTree,
                               String pkFileName, String sk_aaFileName, String ctFileName)  {
        Pairing bp = PairingFactory.getPairing(pairingParametersFileName);

        Properties pkProp = loadPropFromFile(pkFileName);
        Properties sk_aaProp = loadPropFromFile(sk_aaFileName);

        //resolve information from pkFileName
        String gString = pkProp.getProperty("g");
        Element g = bp.getG1().newElementFromBytes(Base64.getDecoder().decode(gString)).getImmutable();
        String hString = pkProp.getProperty("h");
        Element h = bp.getGT().newElementFromBytes(Base64.getDecoder().decode(hString)).getImmutable();
        String egg_alpha_betaString = pkProp.getProperty("egg_alpha_beta");
        Element egg_alpha_beta = bp.getGT().newElementFromBytes(Base64.getDecoder().decode(egg_alpha_betaString)).getImmutable();

        Properties ctProp = new Properties();

        Element s = bp.getZr().newRandomElement().getImmutable();
//        //test
//        Properties skProp = loadPropFromFile(skFileName);
//        String egg_alpha_thetaString = skProp.getProperty("D3");
//        Element egg_alpha_theta = bp.getGT().newElementFromBytes(Base64.getDecoder().decode(egg_alpha_thetaString)).getImmutable();
//        System.out.println("s in the recovered form:"+egg_alpha_theta.powZn(s)+"\n");

        Element omega = bp.getZr().newRandomElement().getImmutable();

        //compute the ciphertext component C3
        Element C3 = message.duplicate().mul(egg_alpha_beta.powZn(omega).getImmutable());

        Element C0 = g.powZn(omega).getImmutable();
        Element C1 = g.powZn(s).getImmutable();
        Element C2 = h.powZn(s).getImmutable();
        Element C1_0 = (C1.mul(C0)).getImmutable();
        Element C2_0 = (C2.mul(h.powZn(omega))).getImmutable();

        ctProp.setProperty("C0", Base64.getEncoder().withoutPadding().encodeToString(C0.toBytes()));
        ctProp.setProperty("C1", Base64.getEncoder().withoutPadding().encodeToString(C1.toBytes()));
        ctProp.setProperty("C2", Base64.getEncoder().withoutPadding().encodeToString(C2.toBytes()));
        ctProp.setProperty("C3", Base64.getEncoder().withoutPadding().encodeToString(C3.toBytes()));
        ctProp.setProperty("C1_0", Base64.getEncoder().withoutPadding().encodeToString(C1_0.toBytes()));
        ctProp.setProperty("C2_0", Base64.getEncoder().withoutPadding().encodeToString(C2_0.toBytes()));

        //set secretshare of root as s
        accessTree[0].secretShare = s;
        //share the secret to all nodes of the accesstree
        nodeShare(accessTree, accessTree[0], bp);

        for (Node node:accessTree) {
            if (node.isLeaf()){

                //resolve information from sk_aaFileName
                //att in the accessTreeList must be a member of the authorityAttList
                String D_aa_String = sk_aaProp.getProperty("D_aa_"+node.att);
                Element D_aa_ = bp.getG1().newElementFromBytes(Base64.getDecoder().decode(D_aa_String)).getImmutable();
                Element C = D_aa_.powZn(node.secretShare);
                ctProp.setProperty("C"+node.att, Base64.getEncoder().withoutPadding().encodeToString(C.toBytes()));
            }
        }
        storePropToFile(ctProp, ctFileName);
    }

    public static Element Decrypt(String pairingParametersFileName, Node[] accessTree, String ctFileName, String sk_aaFileName, String skFileName) {
        Pairing bp = PairingFactory.getPairing(pairingParametersFileName);

        Properties ctProp = loadPropFromFile(ctFileName);
        Properties skProp = loadPropFromFile(skFileName);
        Properties sk_aaProp = loadPropFromFile(sk_aaFileName);

        //resolve information from skprop
        String userAttListString = skProp.getProperty("userAttList");

        //Restore user attribute list int [] type
        int[] userAttList =
                Arrays.stream(userAttListString.substring(1, userAttListString.length()-1).split(",")).map(String::trim).mapToInt(Integer::parseInt).toArray();

        //System.out.println("用户属性列表：" + userAttListString);

        //resolve information from ctprop
        String C1_0String = ctProp.getProperty("C1_0");
        Element C1_0 = bp.getG1().newElementFromBytes(Base64.getDecoder().decode(C1_0String)).getImmutable();
        String C2_0String = ctProp.getProperty("C2_0");
        Element C2_0 = bp.getG1().newElementFromBytes(Base64.getDecoder().decode(C2_0String)).getImmutable();
        String C0String = ctProp.getProperty("C0");
        Element C0 = bp.getG1().newElementFromBytes(Base64.getDecoder().decode(C0String)).getImmutable();
        String C3String = ctProp.getProperty("C3");
        Element C3 = bp.getGT().newElementFromBytes(Base64.getDecoder().decode(C3String)).getImmutable();

        //resolve information from skprop
        String DString = skProp.getProperty("D");
        Element D = bp.getG1().newElementFromBytes(Base64.getDecoder().decode(DString)).getImmutable();
        String D1String = skProp.getProperty("D1");
        Element D1 = bp.getG1().newElementFromBytes(Base64.getDecoder().decode(D1String)).getImmutable();
        String D2String = skProp.getProperty("D2");
        Element D2 = bp.getG1().newElementFromBytes(Base64.getDecoder().decode(D2String)).getImmutable();
//        String D3String = skProp.getProperty("D3");
//        Element D3 = bp.getGT().newElementFromBytes(Base64.getDecoder().decode(D3String)).getImmutable();


        for (Node node : accessTree) {
            if (node.isLeaf()) {
                if (Arrays.stream(userAttList).boxed().toList().contains(node.att)){
                    //only if the node in the userAttList can be found in the accessTree, we have the following computations
                    //that come to the secretShare:
//                    String D_aa_String = sk_aaProp.getProperty("D_aa_"+node.att);
//                    Element D_aa_ = bp.getG1().newElementFromBytes(Base64.getDecoder().decode(D_aa_String)).getImmutable();
                    String CString = ctProp.getProperty("C"+node.att);
                    Element C = bp.getG1().newElementFromBytes(Base64.getDecoder().decode(CString)).getImmutable();
                    String D_usrString = skProp.getProperty("D_usr"+node.att);
                    Element D_usr = bp.getG1().newElementFromBytes(Base64.getDecoder().decode(D_usrString)).getImmutable();

                    node.secretShare = bp.pairing(D_usr,C);
                }
            }
        }
        // nodeRecover to get s
        // s is recovered in the form of egg_alpha_theta_s
        boolean treeOK = nodeRecover(accessTree, accessTree[0], userAttList, bp);
        if (treeOK) {
            //System.out.println("recovered s:"+accessTree[0].secretShare+"\n");
            System.out.println("The access tree is satisfied.");
            Element A = bp.pairing(D1,C1_0).div((accessTree[0].secretShare).mul(bp.pairing(D2,C2_0)));
            Element m=C3.mul(A).div(bp.pairing(D,C0));
            return m;
        }
        else {
            System.out.println("The access tree is not satisfied.");
            return null;
        }
    }

    //q(x)=coef[0] + coef[1]*x^1 + coef[2]*x^2 + coef[d-1]*x^(d-1)
    //hint: the type of exponential is Zr Element
    //randomly select coef[] and make sure that q(0)=s
    public static Element[] randomP(int d, Element s, Pairing bp) {
        Element[] coef = new Element[d];
        coef[0] = s;
        for (int i = 1; i < d; i++){
            coef[i] = bp.getZr().newRandomElement().getImmutable();
        }
        return  coef;
    }
    //to calculate qx(index)
    public static Element qx(Element index, Element[] coef, Pairing bp){
        Element res = coef[0].getImmutable();
        for (int i = 1; i < coef.length; i++){
            //exponential[]
            Element exp = bp.getZr().newElement(i).getImmutable();
            //Note: items that are reused in the loop must be copied out with duplicate
            //getImmutable() can be an option of duplicate()
            res = res.add(coef[i].mul(index.duplicate().powZn(exp)));
        }
        return res;
    }
    //x:target i ∈ S
    public static Element lagrange(int i, int[] S, int x, Pairing bp) {
        Element res = bp.getZr().newOneElement().getImmutable();
        Element iElement = bp.getZr().newElement(i).getImmutable();
        Element xElement = bp.getZr().newElement(x).getImmutable();
        for (int j : S) {
            if (i != j) {
                Element numerator = xElement.sub(bp.getZr().newElement(j));
                Element denominator = iElement.sub(bp.getZr().newElement(j));
                res = res.mul(numerator.div(denominator));
            }
        }
        return res;
    }

    // n is the start point where its secret will be shared
    public static void nodeShare(Node[] nodes, Node n, Pairing bp){
        // only nonleaf nodes should take action
        if (!n.isLeaf()){
            // If it is not a leaf node, a random polynomial is formed, and the constant term of the polynomial
            // is the secret value of the current node (this value will be used for sharing).
            // The polynomial degree is determined by the threshold value gate[0].
            Element[] coef = randomP(n.gate[0], n.secretShare, bp);
            for (int j=0; j<n.children.length; j++ ){
                Node childNode = nodes[n.children[j]];
                childNode.secretShare = qx(bp.getZr().newElement(n.children[j]), coef, bp);
                // Recursion, continue to share the secret of the child node
                nodeShare(nodes, childNode, bp);
            }
        }
    }

    //decide whether node n can be recovered
    public static boolean nodeRecover(Node[] nodes, Node n,  int[] atts, Pairing bp) {
        if (!n.isLeaf()) {
            // For internal nodes, maintain a child node index list for secret recovery.
            List<Integer> validChildrenList = new ArrayList<>();
            int[] validChildren;
            //traverse child nodes of node n
            //the default value of n.valid is false
            for (int j=0; j<n.children.length; j++){
                Node childNode = nodes[n.children[j]];
                if (nodeRecover(nodes, childNode, atts, bp))
                {
                    //System.out.println("The node with index " + n.children[j] + " is satisfied!");
                    //store the index
                    validChildrenList.add(n.children[j]);
                    // If the number of child nodes that meet the conditions has reached the threshold value, the cycle
                    // will jump out and the remaining nodes will not be calculated.
                    if (validChildrenList.size() == n.gate[0]) {
                        n.valid = true;
                        break;
                    }
                }
                else {
                    System.out.println("The node with index " + n.children[j] + " is not satisfied!");
                }
            }
            if (validChildrenList.size() == n.gate[0]){
                validChildren = validChildrenList.stream().mapToInt(i->i).toArray();
                //Note that the Lagrange d-value is made on the [exponential factor] here
                //that is because the recovered secret value is in the form of egg_alpha_theta_qx(0)
                Element secret = bp.getGT().newZeroElement().getImmutable(); //identity element
                for (int i : validChildren) {
                    Element delta = lagrange(i, validChildren, 0, bp);
                    //note that here we do continuously multiply
                    secret = secret.mul(nodes[i].secretShare.duplicate().powZn(delta));
                }
                n.secretShare = secret;
            }
        }
        else {
            //decide whether the attribute value of the leaf node belongs to the attribute list
            if (Arrays.stream(atts).boxed().toList().contains(n.att)){
                n.valid = true;
            }
        }
        return n.valid;
    }

    public static void storePropToFile(Properties prop, String fileName){
        try(FileOutputStream out = new FileOutputStream(fileName)){
            prop.store(out, null);
        }
        catch (IOException e) {
            e.printStackTrace();
            System.out.println(fileName + " save failed!");
            System.exit(-1);
        }
    }

    public static Properties loadPropFromFile(String fileName) {
        Properties prop = new Properties();
        try (FileInputStream in = new FileInputStream(fileName)){
            prop.load(in);
        }
        catch (IOException e){
            e.printStackTrace();
            System.out.println(fileName + " load failed!");
            System.exit(-1);
        }
        return prop;
    }


    public static void main(String[] args) {

        int num=16;
        int[] authorityAttList=new int[num];
        for(int i=0;i<num;i++)
            authorityAttList[i]=i+1;

        /*
        //4 attributes
        int[] userAttList = {1, 2, 3, 4};
        Node[] accessTree = new Node[5];
        accessTree[0] = new Node(new int[]{4,4}, new int[]{1,2,3,4});
        accessTree[1] = new Node(1);
        accessTree[2] = new Node(2);
        accessTree[3] = new Node(3);
        accessTree[4] = new Node(4);

         */

        /*
        //8 attributes:
        int[] userAttList={1,2,3,4,5,6,7,8};
        Node[] accessTree = new Node[14];
        accessTree[0] = new Node(new int[]{2,3}, new int[]{1,2,3});
        accessTree[1] = new Node(new int[]{2,2}, new int[]{4,5});
        accessTree[2] = new Node(new int[]{2,2}, new int[]{6,7});
        accessTree[3] = new Node(new int[]{2,2}, new int[]{8,9});
        accessTree[4] = new Node(1);
        accessTree[5] = new Node(2);
        accessTree[6] = new Node(new int[]{2,2}, new int[]{10,11});
        accessTree[7] = new Node(3);
        accessTree[8] = new Node(4);
        accessTree[9] = new Node(new int[]{2,2},new int[]{12,13});
        accessTree[10] = new Node(5);
        accessTree[11] = new Node(6);
        accessTree[12] = new Node(7);
        accessTree[13] = new Node(8);

         */

        /*
        //12 attributes:
        int[] userAttList={1,2,3,4,5,6,7,8,9,10,11,12};
        Node[] accessTree = new Node[22];
        accessTree[0]= new Node(new int[]{3,4}, new int[]{1,2,3,4});
        accessTree[1]= new Node(new int[]{2,2}, new int[]{5,6});
        accessTree[2]= new Node(new int[]{1,1}, new int[]{7});
        accessTree[3]= new Node(new int[]{2,2}, new int[]{8,9});
        accessTree[4]= new Node(new int[]{2,2}, new int[]{10,11});
        accessTree[6]= new Node(new int[]{2,2}, new int[]{12,13});
        accessTree[8]= new Node(new int[]{2,2}, new int[]{14,15});
        accessTree[9]= new Node(new int[]{2,2}, new int[]{16,17});
        accessTree[11]= new Node(new int[]{2,2}, new int[]{18,19});
        accessTree[16]= new Node(new int[]{2,2}, new int[]{20,21});
        accessTree[5] = new Node(1);
        accessTree[7] = new Node(4);
        accessTree[10] = new Node(10);
        accessTree[12] = new Node(2);
        accessTree[13] = new Node(3);
        accessTree[14] = new Node(5);
        accessTree[15] = new Node(6);
        accessTree[17] = new Node(9);
        accessTree[18] = new Node(11);
        accessTree[19] = new Node(12);
        accessTree[20] = new Node(7);
        accessTree[21] = new Node(8);

         */

        //16 attributes:
        int[] userAttList={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
        Node[] accessTree = new Node[28];
        accessTree[0]= new Node(new int[]{4,5}, new int[]{1,2,3,4,5});
        accessTree[1]= new Node(new int[]{2,2}, new int[]{6,7});
        accessTree[2]= new Node(new int[]{2,2}, new int[]{8,9});
        accessTree[3]= new Node(new int[]{2,2}, new int[]{10,11});
        accessTree[4]= new Node(new int[]{2,2}, new int[]{12,13});
        accessTree[5]= new Node(new int[]{1,1}, new int[]{14});
        accessTree[7]= new Node(new int[]{2,3}, new int[]{15,16,17});
        accessTree[8]= new Node(new int[]{1,1}, new int[]{18});
        accessTree[9]= new Node(new int[]{2,2}, new int[]{19,20});
        accessTree[12]= new Node(new int[]{2,3}, new int[]{21,22,23});
        accessTree[14]= new Node(new int[]{2,2}, new int[]{24,25});
        accessTree[20]= new Node(new int[]{2,2}, new int[]{26,27});
        accessTree[6] = new Node(1);
        accessTree[10] = new Node(9);
        accessTree[11] = new Node(10);
        accessTree[13] = new Node(14);
        accessTree[15] = new Node(2);
        accessTree[16] = new Node(3);
        accessTree[17] = new Node(4);
        accessTree[18] = new Node(5);
        accessTree[19] = new Node(6);
        accessTree[21] = new Node(11);
        accessTree[22] = new Node(12);
        accessTree[23] = new Node(13);
        accessTree[24] = new Node(15);
        accessTree[25] = new Node(16);
        accessTree[26] = new Node(7);
        accessTree[27] = new Node(8);


        //automatically go to root dir to find files
        String dir = "data/";
        String pairingParametersFileName = "a.properties";
        String pkFileName = dir + "pk.properties";
        String mkFileName = dir + "mk.properties";
        String skFileName = dir + "sk.properties";
        String ctFileName = dir + "ct.properties";
        String sk_aaFileName = dir + "sk_aa.properties";

        long starttime=System.currentTimeMillis();
        setup(pairingParametersFileName, pkFileName, mkFileName);
        aa_keygen(pairingParametersFileName,authorityAttList, pkFileName, mkFileName, sk_aaFileName);
        long endtime=System.currentTimeMillis();
        System.out.println("初始化用时：" + (endtime-starttime)+"ms\n");

        starttime=System.currentTimeMillis();
        keygen(pairingParametersFileName, userAttList, pkFileName, mkFileName, sk_aaFileName, skFileName);
        endtime=System.currentTimeMillis();
        System.out.println("密钥生成用时：" + (endtime-starttime)+"ms\n");

        Element message = PairingFactory.getPairing(pairingParametersFileName).getGT().newRandomElement().getImmutable();
        //System.out.println("明文消息:" + message);

        starttime=System.currentTimeMillis();
        encrypt(pairingParametersFileName, message, accessTree, pkFileName, sk_aaFileName, ctFileName);
        endtime=System.currentTimeMillis();
        System.out.println("加密用时：" + (endtime-starttime)+"ms\n");

        starttime=System.currentTimeMillis();
        Element res = Decrypt(pairingParametersFileName, accessTree, ctFileName, sk_aaFileName, skFileName);
        endtime=System.currentTimeMillis();
        System.out.println("解密用时：" + (endtime-starttime)+"ms\n");

        //System.out.println("解密结果:" + res);

        if (message.isEqual(res)) {
            System.out.println("成功解密！");
        }
    }

}

