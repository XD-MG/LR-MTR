function [ rate1, rate2, rate3, rate4,rate5] = RATE_Classfier(predict_label, test,true_testlabel,NUM)
     n1=NUM(1);n2=NUM(2);n3=NUM(3);n4=NUM(4);n5=NUM(5);   
testlabel=predict_label;
     num=0; num1=0; num2=0; num3=0; num4=0;
     
     for i  = 1:size(test,2)
     AA = testlabel(i);
     BB = true_testlabel(i);
     
    if AA == BB && BB == 1
    num = num+1;
    end
    
    if AA == BB && BB == 2
    num1 = num1+1;
    end
    
    if AA == BB && BB == 3
    num2 = num2+1;
    end
    
    if AA == BB && BB == 4
    num3 = num3+1;
    end
    
    if AA == BB && BB == 5
    num4 = num4+1;
    end
     end
   rate1=(num/n1)*100;
    rate2=(num1/n2)*100;
    rate3=(num2/n3)*100;
    rate4=(num3/n4)*100;
    rate5=(num4/n5)*100;
end