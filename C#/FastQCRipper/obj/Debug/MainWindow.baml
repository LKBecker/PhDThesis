   M S B A M L   `   `   `  ����G  CFastQCRipper, Version=1.0.0.0, Culture=neutral, PublicKeyToken=null+clr-namespace:FastQCRipperFastQCRipper  5   	       FastQCRipper.MainWindow   5   	   W SPresentationCore, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35R NWindowsBase, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35\ XPresentationFramework, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35D 9http://schemas.microsoft.com/winfx/2006/xaml/presentation    5   	   R NSystem.Xaml, Version=4.0.0.0, Culture=neutral, PublicKeyToken=b77a5c561934e089:x,http://schemas.microsoft.com/winfx/2006/xaml     5   	   @d2http://schemas.microsoft.com/expression/blend/2008     5   	   Jmc;http://schemas.openxmlformats.org/markup-compatibility/2006     5   	   ,localclr-namespace:FastQCRipper     5   	     � Title$  Fast QC Ripper��5   	    � SizeToContent$ WidthAndHeight=�6     � 
ResizeMode$ CanResizeWithGrip=�6?   $��Center=�6^   $��Center=�6{   .��5	      � -   5         �Name$ MainGrid��5	      $	��521��6   $	��105��6&    J� RenderTransformOrigin$ 0.5,0.5&�63   $	��Top=�6S   .+�5
      �� -   5         ��Name$! ButtonSelectSourceDirectory��5
      $��	Select...��60   $
��Left=�6D   $��	436,9,0,0q�6_   $	��Top=�6r   $��75��6�   6�   �� 5      -   5         ��Name$ TextBoxSourceDirectory��5      $
��Left=�6,   $��23��6G   $��	10,10,0,0q�6S    �� TextWrapping$
 Wrap=�6f   $��(No directory selected)��6z   $	��Top=�6�   $	��400��6�   4  .System.Windows.Automation.AutomationProperties   Name$ TextBox.SourceDirectory��6�   	   IsRowHeader	 . 
 �� 
IsReadOnly
 .  �� MaxLines$ 1��6+   �� MinLines$ 1��68  
���   �6E  �� 5      -   5        $&  ButtonSelectDestinationDirectory��5      $��	Select...��65   $
��Left=�6I   $��
436,37,0,0q�6d   $	��Top=�6x   $��75��6�   6�   �� 5      -   5        $! TextBoxDestinationDirectory��5      $
��Left=�61   $��23��6L   $��	10,38,0,0q�6X   $
 Wrap=�6k   $��(No directory selected)��6   $	��Top=�6�   $	��400��6�   $ 1��6�   
 . 
���   �6�   �� 5      -   5        $ ButtomCompute��5      $��Extract��6"   $
��Left=�64   $��
436,75,0,0q�6O   $	��Top=�6c   $��75��6{   6�   � 5      -   5         �Name$ ProgressBarZIPRip��5      $��Center=�6+   $��19��6H   $��0q�6T   $	��400��6_   $��Center=�6k   $y�Hidden=�6�   �� 5      -   5         ��Name$ LabelProgress��5   
   $'��!Processing File X of Y (ZZZ%) ...��6!   $��Center=�6M   $��0,65,0,0q�6j   $	��Top=�6|   $y�Hidden=�6�   9� 5      $��Options��6   $��61��6   $��10,75,111,-31q�6*   $��Bottom=�6A    9� 
IsExpanded .  .��5      �� $
��Left=�6   $��51��6+   $��10,0,0,0q�67   $	��380��6I    �� Orientation$ 
Horizontal=�6U   .+�5      �� $��File search pattern��6   $��23��6*   $��0q�66   $	��110��6A   $
��Left=�6M   $��Center=�6h   �� 5      $
��Left=�6   $��23��6)   $
 Wrap=�65   $��TextBox��6I   $��Center=�6X   $	��240��6s   $��
10,14,0,14q�6   5      5      5      5      